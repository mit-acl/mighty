/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "dgp/dgp_manager.hpp"

/// The type of map data Tmap is defined as a 1D array
using Tmap = std::vector<char>;
using namespace mighty;
using namespace termcolor;

typedef timer::Timer MyTimer;

DGPManager::DGPManager() {}

void DGPManager::setParameters(const parameters &par)
{

    // Get the parameter
    par_ = par;

    // Set the parameters
    res_ = par.res;
    drone_radius_ = par.drone_radius;
    max_dist_vertexes_ = par.max_dist_vertexes;
    use_shrinked_box_ = par.use_shrinked_box;
    shrinked_box_size_ = par.shrinked_box_size;
    // local_box_size_ is std::vector<float> but par.local_box_size is std::vector<double>
    local_box_size_ = {static_cast<float>(par.local_box_size[0]), static_cast<float>(par.local_box_size[1]), static_cast<float>(par.local_box_size[2])};

    // shared pointer to the map util for actual planning
    map_util_ = std::make_shared<mighty::VoxelMapUtil>(par.factor_dgp * par.res, par.x_min, par.x_max, par.y_min, par.y_max, par.z_min, par.z_max, par.inflation_dgp);
}

void DGPManager::cleanUpPath(vec_Vecf<3> &path)
{
    planner_ptr_->cleanUpPath(path);
}

void DGPManager::setupDGPPlanner(const std::string &global_planner, bool global_planner_verbose, double res, double v_max, double a_max, double j_max, int dgp_timeout_duration_ms, double w_unknown, double w_align, double decay_len_cells, double w_side, int los_cells, double min_len, double min_turn)
{

    // Get the parameters
    v_max_ = v_max;
    v_max_3d_ = Eigen::Vector3d(v_max, v_max, v_max);
    a_max_ = a_max;
    a_max_3d_ = Eigen::Vector3d(a_max, a_max, a_max);
    j_max_ = j_max;
    j_max_3d_ = Eigen::Vector3d(j_max, j_max, j_max);

    // Create the DGP planner
    planner_ptr_ = std::unique_ptr<DGPPlanner>(new DGPPlanner(global_planner, global_planner_verbose, v_max, a_max, j_max, dgp_timeout_duration_ms, w_unknown, w_align, decay_len_cells, w_side, los_cells, min_len, min_turn));

    // Create the map_util_for_planning
    // This is the beginning of the planning, so we fetch the map_util_ and don't update it for the entire planning process (updating while planning makes the planner slower)
    mtx_map_util_.lock();
    map_util_for_planning_ = std::make_shared<mighty::VoxelMapUtil>(*map_util_);
    mtx_map_util_.unlock();
}

void DGPManager::updateVmax(double v_max)
{
    // Update the maximum velocity
    v_max_ = v_max;
    v_max_3d_ = Eigen::Vector3d(v_max, v_max, v_max);
    planner_ptr_->updateVmax(v_max);
}

void DGPManager::freeStart(Vec3f &start, double factor)
{
    // Set start free
    Veci<3> start_int = map_util_for_planning_->floatToInt(start);
    map_util_for_planning_->setFreeVoxelAndSurroundings(start_int, factor * res_);
    // map_util_for_planning_->setFree(start_int);
}

void DGPManager::freeGoal(Vec3f &goal, double factor)
{
    // Set goal free
    Veci<3> goal_int = map_util_for_planning_->floatToInt(goal);
    map_util_for_planning_->setFreeVoxelAndSurroundings(goal_int, factor * res_);
    // map_util_for_planning_->setFree(goal_int);
}

bool DGPManager::checkIfPointOccupied(const Vec3f &point)
{
    // Check if the point is free
    Veci<3> point_int = map_util_for_planning_->floatToInt(point);

    return map_util_for_planning_->isOccupied(point_int);
}

// Sample along [p0, p1] at a safe step to ensure we don't skip thin obstacles.
// Uses the occupancy from the (already inflated) planning map.
inline bool isSegmentFree(const mighty::VoxelMapUtil &map,
                          const Vec3f &p0,
                          const Vec3f &p1,
                          const double sample_step)
{
    const Vec3f d = p1 - p0;
    const double L = d.norm();
    if (L <= 1e-6)
        return true;

    const Vec3f dir = d / L;
    // March from start to end, including the goal voxel.
    for (double s = 0.0; s <= L; s += sample_step)
    {
        const Vec3f q = p0 + s * dir;
        const Veci<3> qi = map.floatToInt(q);
        // if (!map.isFree(qi))
        if (map.isOccupied(qi))
            return false;
    }
    // Ensure exact endpoint is also checked (if s stepped past)
    const Veci<3> q_end = map.floatToInt(p1);
    // if (!map.isFree(q_end))
    if (map.isOccupied(q_end))
        return false;

    return true;
}

// Greedily collapse a path into maximal collision-free segments.
// This mirrors the "generate a long segment if it’s collision free" behavior.
inline void collapseIntoLongSegments(const mighty::VoxelMapUtil &map,
                                     double res,
                                     vec_Vecf<3> &path_inout,
                                     double sample_step = -1.0)
{
    if (path_inout.size() <= 2)
        return;

    const double step = (sample_step > 0.0) ? sample_step : 0.5 * res;

    vec_Vecf<3> simplified;
    simplified.reserve(path_inout.size());
    simplified.push_back(path_inout.front()); // keep start

    size_t anchor = 0;     // current segment start index
    size_t j = anchor + 1; // candidate end

    while (j < path_inout.size())
    {
        size_t last_good = anchor + 1;

        // Extend j as far as LoS holds
        while (j < path_inout.size() &&
               isSegmentFree(map, path_inout[anchor], path_inout[j], step))
        {
            last_good = j;
            ++j;
        }

        // Commit the farthest valid endpoint
        simplified.push_back(path_inout[last_good]);

        // Start next segment from there
        anchor = last_good;
        j = anchor + 1;
    }

    // Make sure we end exactly at the original goal
    simplified.back() = path_inout.back();

    path_inout.swap(simplified);
}

bool DGPManager::solveDGP(const Vec3f &start_sent, const Vec3f &start_vel, const Vec3f &goal_sent, double &final_g, double weight, double current_time, vec_Vecf<3> &path)
{
    // Set start and goal
    Eigen::Vector3d start(start_sent(0), start_sent(1), start_sent(2));
    Eigen::Vector3d goal(goal_sent(0), goal_sent(1), goal_sent(2));

    //
    // Run DGP
    //

    // Set collision checking function
    planner_ptr_->setMapUtil(map_util_for_planning_);

    // DGP Plan
    bool result = false;

    // Attempt to plan
    result = planner_ptr_->plan(start, start_vel, goal, final_g, current_time, weight);

    // If there is a solution
    if (result)
    {
        path = planner_ptr_->getPath();
    }
    else
    {
        // std::cout << bold << red << "DGP didn't find a solution from " << start.transpose() << " to " << goal.transpose() << reset << std::endl;
    }

    // Clean up path
    planner_ptr_->cleanUpPath(path);

    // // Add more vertices if necessary
    mighty_utils::createMoreVertexes(path, max_dist_vertexes_);

    return result;
}

bool DGPManager::checkIfPathInFree(const vec_Vecf<3> &path, vec_Vecf<3> &free_path)
{
    // Initialize result
    free_path.clear();
    free_path.push_back(path[0]);

    // Plan only in free space if required
    for (size_t i = 1; i < path.size(); i++)
    {
        Veci<3> path_int = map_util_for_planning_->floatToInt(path[i]);
        if (map_util_for_planning_->isFree(path_int))
        {
            free_path.push_back(path[i]);
        }
        else
        {
            break;
        }
    }

    if (free_path.size() <= 1)
        return false;

    return true;
}

void DGPManager::pushPathIntoFreeSpace(const vec_Vecf<3> &path, vec_Vecf<3> &free_path)
{

    // Initialize result
    free_path.clear();
    free_path.push_back(path[0]);

    // Plan only in free space if required
    for (size_t i = 1; i < path.size(); i++)
    {
        Vec3f free_point;
        map_util_for_planning_->findClosestFreePoint(path[i], free_point);
        free_path.push_back(free_point);
    }
}

bool DGPManager::checkIfPointFree(const Vec3f &point) const
{
    // Check if the point is free
    Veci<3> point_int = map_util_for_planning_->floatToInt(point);
    return map_util_for_planning_->isFree(point_int);
}

bool DGPManager::checkIfPointHasNonFreeNeighbour(const Vec3f &point) const
{
    // Check if the point has an occupied neighbour
    Veci<3> point_int = map_util_for_planning_->floatToInt(point);
    return map_util_for_planning_->checkIfPointHasNonFreeNeighbour(point_int);
}

void DGPManager::getOccupiedCells(vec_Vecf<3> &occupied_cells)
{
    // Get the occupied cells
    mtx_map_util_.lock();
    occupied_cells = map_util_->getOccupiedCloud();
    mtx_map_util_.unlock();
}

void DGPManager::getFreeCells(vec_Vecf<3> &free_cells)
{
    // Get the free cells
    mtx_map_util_.lock();
    free_cells = map_util_->getFreeCloud();
    mtx_map_util_.unlock();
}

void DGPManager::getComputationTime(double &global_planning_time, double &dgp_static_jps_time, double &dgp_check_path_time, double &dgp_dynamic_astar_time, double &dgp_recover_path_time)
{
    // Get the computation time
    global_planning_time = planner_ptr_->getInitialGuessPlanningTime();
    dgp_static_jps_time = planner_ptr_->getStaticJPSPlanningTime();
    dgp_check_path_time = planner_ptr_->getCheckPathTime();
    dgp_dynamic_astar_time = planner_ptr_->getDynamicAstarTime();
    dgp_recover_path_time = planner_ptr_->getRecoverPathTime();
}

void DGPManager::getVecOccupied(vec_Vec3f &vec_o)
{
    mtx_vec_o_.lock();
    vec_o = vec_o_;
    mtx_vec_o_.unlock();
}

void DGPManager::updateVecOccupied(const vec_Vec3f &vec_o)
{
    mtx_vec_o_.lock();
    vec_o_ = vec_o;
    mtx_vec_o_.unlock();
}

void DGPManager::getVecUnknownOccupied(vec_Vec3f &vec_uo)
{
    mtx_vec_uo_.lock();
    vec_uo = vec_uo_;
    mtx_vec_uo_.unlock();
}

void DGPManager::updateVecUnknownOccupied(const vec_Vec3f &vec_uo)
{
    mtx_vec_uo_.lock();
    vec_uo_ = vec_uo;
    mtx_vec_uo_.unlock();
}

void DGPManager::insertVecOccupiedToVecUnknownOccupied()
{
    mtx_vec_uo_.lock();
    mtx_vec_o_.lock();
    vec_uo_.insert(vec_uo_.end(), vec_o_.begin(), vec_o_.end());
    mtx_vec_o_.unlock();
    mtx_vec_uo_.unlock();
}

bool DGPManager::cvxEllipsoidDecomp(const state &A, const vec_Vecf<3> &path,
                                    std::vector<LinearConstraint3D> &l_constraints,
                                    vec_E<Polyhedron<3>> &poly_out,
                                    bool use_for_safe_path)
{

    // Initialize result.
    bool result = true;

    // Get unknown occupied cells.
    if (use_for_safe_path)
    {
        // If we are using the convex decomposition for safe paths, we include both unknown and occupied cells.
        vec_Vecf<3> vec_uo;
        getVecUnknownOccupied(vec_uo);
        ellip_decomp_util_.set_obs(vec_uo);
    }
    else
    {
        // Otherwise, we only include occupied cells.
        vec_Vecf<3> vec_o;
        getVecOccupied(vec_o);
        ellip_decomp_util_.set_obs(vec_o);
    }

    // Set the local bounding box and z constraints.
    ellip_decomp_util_.set_local_bbox(Vec3f(local_box_size_[0], local_box_size_[1], local_box_size_[2]));
    ellip_decomp_util_.set_z_min_and_max(par_.z_min, par_.z_max); // buffer for the drone size

    // Set the inflate distance.
    ellip_decomp_util_.set_inflate_distance(drone_radius_);

    // Find convex polyhedra.
    ellip_decomp_util_.dilate(path, result);

    if (!result)
        return false;
    // Optionally shrink polyhedra.
    if (use_shrinked_box_)
        ellip_decomp_util_.shrink_polyhedrons(shrinked_box_size_);

    // Get the polyhedra.
    auto polys = ellip_decomp_util_.get_polyhedrons();

    // Preallocate the constraints vector.
    size_t numConstraints = (path.size() > 0) ? (path.size() - 1) : 0;
    l_constraints.clear();
    l_constraints.resize(numConstraints);

    // Flag to record if any thread finds an error.
    bool errorFound = false;

// Parallelize the constraint computation loop.
#pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(numConstraints); i++)
    {

        // Compute the midpoint between consecutive path points.
        auto pt_inside = (path[i] + path[i + 1]) / 2.0;
        LinearConstraint3D cs(pt_inside, polys[i].hyperplanes(), polys[i]);

        // If either matrix A_ or vector b_ contains NaN, mark an error.
        if (cs.A_.hasNaN() || cs.b_.hasNaN())
        {
#pragma omp atomic write
            errorFound = true;
        }
        else
        {
            l_constraints[i] = cs;
        }
    }

    // If an error was detected, report and exit.
    if (errorFound)
    {
        std::cout << "A_ or b_ has NaN" << std::endl;
        return false;
    }

    // Return the computed polyhedra.
    poly_out = std::move(polys);

    return true;
}

void DGPManager::updateMap(double wdx, double wdy, double wdz, const Vec3f &center_map, const pcl::PointCloud<pcl::PointXYZ>::ConstPtr &pclptr)
{

    // Get the current time to see the computation time for readmap
    auto start_time = std::chrono::high_resolution_clock::now();

    mtx_map_util_.lock();
    map_util_->readMap(pclptr, (int)wdx / res_, (int)wdy / res_, (int)wdz / res_, center_map, par_.z_min, par_.z_max, par_.inflation_dgp); // Map read
    mtx_map_util_.unlock();

    // Get the elapsed time for reading the map
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    // std::cout << "Map read time: " << elapsed_time << " ms" << std::endl;

    if (!map_initialized_)
    {
        map_initialized_ = true;
    }
}

bool DGPManager::isMapInitialized() const
{
    return map_initialized_;
}

void DGPManager::findClosestFreePoint(const Vec3f &point, Vec3f &closest_free_point)
{
    mtx_map_util_.lock();
    map_util_->findClosestFreePoint(point, closest_free_point);
    mtx_map_util_.unlock();
}

int DGPManager::countUnknownCells() const
{
    return map_util_for_planning_->countUnknownCells();
}

int DGPManager::getTotalNumCells() const
{
    return map_util_for_planning_->getTotalNumCells();
}