/*
    MIT License

    Copyright (c) 2021 Zhepei Wang (wangzhepei@live.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef SFC_GEN_HPP
#define SFC_GEN_HPP

#include "geo_utils.hpp"
#include "firi.hpp"

#include <ompl/util/Console.h>
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/geometric/planners/rrt/InformedRRTstar.h>
#include <ompl/base/objectives/PathLengthOptimizationObjective.h>
#include <ompl/base/DiscreteMotionValidator.h>

#include <deque>
#include <memory>
#include <Eigen/Eigen>

namespace sfc_gen
{

    template <typename Map>
    inline double planPath(const Eigen::Vector3d &s,
                           const Eigen::Vector3d &g,
                           const Eigen::Vector3d &lb,
                           const Eigen::Vector3d &hb,
                           Map *mapPtr,
                           const double &timeout,
                           std::vector<Eigen::Vector3d> &p)
    {
        auto space(std::make_shared<ompl::base::RealVectorStateSpace>(3));

        ompl::base::RealVectorBounds bounds(3);
        bounds.setLow(0, 0.0);
        bounds.setHigh(0, hb(0) - lb(0));
        bounds.setLow(1, 0.0);
        bounds.setHigh(1, hb(1) - lb(1));
        bounds.setLow(2, 0.0);
        bounds.setHigh(2, hb(2) - lb(2));
        space->setBounds(bounds);

        auto si(std::make_shared<ompl::base::SpaceInformation>(space));

        si->setStateValidityChecker(
            [&](const ompl::base::State *state)
            {
                const auto *pos = state->as<ompl::base::RealVectorStateSpace::StateType>();
                const Eigen::Vector3d position(lb(0) + (*pos)[0],
                                               lb(1) + (*pos)[1],
                                               lb(2) + (*pos)[2]);
                return mapPtr->query(position) == 0;
            });
        si->setup();

        ompl::msg::setLogLevel(ompl::msg::LOG_NONE);

        ompl::base::ScopedState<> start(space), goal(space);
        start[0] = s(0) - lb(0);
        start[1] = s(1) - lb(1);
        start[2] = s(2) - lb(2);
        goal[0] = g(0) - lb(0);
        goal[1] = g(1) - lb(1);
        goal[2] = g(2) - lb(2);

        auto pdef(std::make_shared<ompl::base::ProblemDefinition>(si));
        pdef->setStartAndGoalStates(start, goal);
        pdef->setOptimizationObjective(std::make_shared<ompl::base::PathLengthOptimizationObjective>(si));
        auto planner(std::make_shared<ompl::geometric::InformedRRTstar>(si));
        planner->setProblemDefinition(pdef);
        planner->setup();

        ompl::base::PlannerStatus solved;
        solved = planner->ompl::base::Planner::solve(timeout);

        double cost = INFINITY;
        if (solved)
        {
            p.clear();
            const ompl::geometric::PathGeometric path_ =
                ompl::geometric::PathGeometric(
                    dynamic_cast<const ompl::geometric::PathGeometric &>(*pdef->getSolutionPath()));
            for (size_t i = 0; i < path_.getStateCount(); i++)
            {
                const auto state = path_.getState(i)->as<ompl::base::RealVectorStateSpace::StateType>()->values;
                p.emplace_back(lb(0) + state[0], lb(1) + state[1], lb(2) + state[2]);
            }
            cost = pdef->getSolutionPath()->cost(pdef->getOptimizationObjective()).value();
        }

        return cost;
    }

    inline void convexCover(const std::vector<Eigen::Vector3d> &path,
                            const std::vector<Eigen::Vector3d> &points,
                            const Eigen::Vector3d &lowCorner,
                            const Eigen::Vector3d &highCorner,
                            const double &progress,
                            const double &range,
                            std::vector<Eigen::MatrixX4d> &hpolys,
                            const double eps = 1.0e-6)
    {
        hpolys.clear();
        const int n = path.size();
        if (n < 2)
            return;

        // --------------------------------------------------------------------
        // Optional global safety buffer: shrink all final planes by this amount.
        // Set to (robot_radius + safety_margin) if you want geometry-only inflation.
        // Leave at 0.0 to preserve previous behavior.
        constexpr double CLEARANCE = 0.0;
        // --------------------------------------------------------------------

        auto shrinkPoly = [](Eigen::MatrixX4d &hp, double clearance)
        {
            if (clearance <= 0.0)
                return;
            for (int r = 0; r < hp.rows(); ++r)
            {
                const double nn = hp(r, 0) * hp(r, 0) + hp(r, 1) * hp(r, 1) + hp(r, 2) * hp(r, 2);
                if (nn > 1e-16)
                    hp(r, 3) += clearance * std::sqrt(nn);
            }
        };

        Eigen::Matrix<double, 6, 4> bd = Eigen::Matrix<double, 6, 4>::Zero();
        bd(0, 0) = 1.0;
        bd(1, 0) = -1.0;
        bd(2, 1) = 1.0;
        bd(3, 1) = -1.0;
        bd(4, 2) = 1.0;
        bd(5, 2) = -1.0;

        Eigen::MatrixX4d hp, gap;
        Eigen::Vector3d a, b = path[0];

        std::vector<Eigen::Vector3d> valid_pc;
        valid_pc.reserve(points.size());
        std::vector<Eigen::Vector3d> bs; // kept for parity with the original
        bs.reserve(static_cast<size_t>(std::ceil((n - 1) * 1.0)));

        for (int i = 1; i < n; /* i is advanced inside */)
        {
            a = b;
            if ((a - path[i]).norm() > progress)
            {
                b = (path[i] - a).normalized() * progress + a;
            }
            else
            {
                b = path[i];
                ++i;
            }
            bs.emplace_back(b);

            // Build local bounding slab
            bd(0, 3) = -std::min(std::max(a(0), b(0)) + range, highCorner(0)); //  x ≤ ...
            bd(1, 3) = +std::max(std::min(a(0), b(0)) - range, lowCorner(0));  //  x ≥ ...
            bd(2, 3) = -std::min(std::max(a(1), b(1)) + range, highCorner(1)); //  y ≤ ...
            bd(3, 3) = +std::max(std::min(a(1), b(1)) - range, lowCorner(1));  //  y ≥ ...
            bd(4, 3) = -std::min(std::max(a(2), b(2)) + range, highCorner(2)); //  z ≤ ...
            bd(5, 3) = +std::max(std::min(a(2), b(2)) - range, lowCorner(2));  //  z ≥ ...

            // Filter obstacle points inside the slab — use robust boundary inclusion
            valid_pc.clear();
            for (const Eigen::Vector3d &p : points)
            {
                // CHANGED: <= eps (was < 0.0). This keeps boundary points.
                if ((bd.leftCols<3>() * p + bd.rightCols<1>()).maxCoeff() <= eps)
                {
                    valid_pc.emplace_back(p);
                }
            }

            // Build obstacle matrix safely (avoid Map on empty vector)
            Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor> pc_mat;
            pc_mat.resize(3, static_cast<int>(valid_pc.size()));
            for (int k = 0; k < static_cast<int>(valid_pc.size()); ++k)
            {
                pc_mat.col(k) = valid_pc[static_cast<size_t>(k)];
            }

            bool used_fallback = false;

            if (!valid_pc.empty())
            {
                // Normal FIRI inflation against local obstacles
                firi::firi(bd, pc_mat, a, b, hp);
            }
            else
            {
                // CHANGED: Robust fallback when the slab has zero obstacle points.
                // Use the slab itself and shrink by CLEARANCE so it stays conservative.
                hp = bd;
                shrinkPoly(hp, CLEARANCE);
                used_fallback = true;
            }

            // CHANGED: Optional global shrink to add a uniform safety buffer.
            shrinkPoly(hp, CLEARANCE);

            // Insert a tiny "gap" poly if both current & previous polys contain 'a'
            if (!hpolys.empty())
            {
                const Eigen::Vector4d ah(a(0), a(1), a(2), 1.0);
                const int cnt = ((hp * ah).array() > -eps).cast<int>().sum() + ((hpolys.back() * ah).array() > -eps).cast<int>().sum();

                if (cnt >= 3)
                {
                    if (!used_fallback)
                    {
                        // Normal connector via FIRI
                        firi::firi(bd, pc_mat, a, a, gap, 1);
                        shrinkPoly(gap, CLEARANCE);
                        hpolys.emplace_back(gap);
                    }
                    else
                    {
                        // CHANGED: If no obstacles, make a conservative slab connector
                        gap = bd;
                        shrinkPoly(gap, CLEARANCE);
                        hpolys.emplace_back(gap);
                    }
                }
            }

            hpolys.emplace_back(hp);
        }
    }

    inline void shortCut(std::vector<Eigen::MatrixX4d> &hpolys)
    {
        // Adaptive, robust post-processing of SFC list.
        // Changes vs. original:
        //  - Adaptive overlap tolerance based on corridor scale (min axis width), clamped in [5mm, 5cm]
        //  - Safer handling of degenerate/empty inputs
        //  - Avoid duplicate consecutive entries in the output

        auto axisExtent = [](const Eigen::MatrixX4d &hp, int axis) -> double
        {
            // Try to read axis-aligned faces that come from the initial bounding slab.
            // We look for normals exactly equal to [±1,0,0], [0,±1,0], [0,0,±1].
            bool hasPos = false, hasNeg = false;
            double pos = 0.0, neg = 0.0;
            const int R = static_cast<int>(hp.rows());
            for (int r = 0; r < R; ++r)
            {
                const double nx = hp(r, 0), ny = hp(r, 1), nz = hp(r, 2), d = hp(r, 3);
                if (axis == 0)
                {
                    if (std::abs(nx - 1.0) < 1e-12 && std::abs(ny) < 1e-12 && std::abs(nz) < 1e-12)
                    {
                        pos = -d;
                        hasPos = true;
                    } // x <= -d
                    if (std::abs(nx + 1.0) < 1e-12 && std::abs(ny) < 1e-12 && std::abs(nz) < 1e-12)
                    {
                        neg = d;
                        hasNeg = true;
                    } // x >=  d
                }
                else if (axis == 1)
                {
                    if (std::abs(ny - 1.0) < 1e-12 && std::abs(nx) < 1e-12 && std::abs(nz) < 1e-12)
                    {
                        pos = -d;
                        hasPos = true;
                    } // y <= -d
                    if (std::abs(ny + 1.0) < 1e-12 && std::abs(nx) < 1e-12 && std::abs(nz) < 1e-12)
                    {
                        neg = d;
                        hasNeg = true;
                    } // y >=  d
                }
                else // axis == 2
                {
                    if (std::abs(nz - 1.0) < 1e-12 && std::abs(nx) < 1e-12 && std::abs(ny) < 1e-12)
                    {
                        pos = -d;
                        hasPos = true;
                    } // z <= -d
                    if (std::abs(nz + 1.0) < 1e-12 && std::abs(nx) < 1e-12 && std::abs(ny) < 1e-12)
                    {
                        neg = d;
                        hasNeg = true;
                    } // z >=  d
                }
            }
            if (hasPos && hasNeg)
                return pos - neg; // width along this axis
            return -1.0;          // invalid
        };

        auto minAxisWidth = [&](const Eigen::MatrixX4d &A, const Eigen::MatrixX4d &B) -> double
        {
            double res = 1e9; // large
            bool found = false;
            for (int ax = 0; ax < 3; ++ax)
            {
                double ea = axisExtent(A, ax);
                double eb = axisExtent(B, ax);
                double e = 1e9;
                bool ok = false;
                if (ea > 0.0)
                {
                    e = std::min(e, ea);
                    ok = true;
                }
                if (eb > 0.0)
                {
                    e = std::min(e, eb);
                    ok = true;
                }
                if (ok)
                {
                    res = std::min(res, e);
                    found = true;
                }
            }
            return found ? res : -1.0;
        };

        auto overlapTol = [&](const Eigen::MatrixX4d &A, const Eigen::MatrixX4d &B) -> double
        {
            // 10% of the corridor's narrowest axis width, clamped to [5mm, 5cm]
            double e = minAxisWidth(A, B);
            double tol = 0.01; // default 1 cm
            if (e > 1e-6)
            {
                tol = 0.1 * e;
                if (tol < 0.005)
                    tol = 0.005;
                else if (tol > 0.05)
                    tol = 0.05;
            }
            return tol;
        };

        if (hpolys.empty())
            return;

        std::vector<Eigen::MatrixX4d> htemp = hpolys;
        if (htemp.size() == 1)
        {
            htemp.insert(htemp.begin(), htemp.front());
        }
        hpolys.clear();

        const int M = static_cast<int>(htemp.size());
        std::deque<int> indices;
        indices.push_front(M - 1);

        for (int i = M - 1; i >= 0; i--)
        {
            for (int j = 0; j < i; ++j)
            {
                bool overlap = (j < i - 1)
                                   ? geo_utils::overlap(htemp[i], htemp[j], overlapTol(htemp[i], htemp[j]))
                                   : true; // always accept immediate neighbor
                if (overlap)
                {
                    indices.push_front(j);
                    i = j + 1; // jump (the outer for will do i-- next)
                    break;
                }
            }
        }

        // Emit, skipping consecutive duplicates
        int last = -1;
        hpolys.reserve(indices.size());
        for (const int idx : indices)
        {
            if (idx == last)
                continue;
            hpolys.emplace_back(htemp[idx]);
            last = idx;
        }
    }

}

#endif
