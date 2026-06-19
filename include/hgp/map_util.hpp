/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

/**
 * @file map_util.h
 * @brief MapUtil classes
 */
#pragma once

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <mutex>

#include <pcl/kdtree/kdtree_flann.h>

#include <mighty/mighty_type.hpp>

#include "hgp/data_type.hpp"

#include "timer.hpp"
#include <mighty/esdf_grid_2d.hpp>
#include <mighty/occ_grid_2d.hpp>

namespace mighty {

// The type of map data Tmap is defined as a 1D array
// Using int8_t saves 75% memory compared to int (only need 3 values: -1, 0, 100)
using Tmap = std::vector<int8_t>;
typedef timer::Timer MyTimer;

/**
 * @brief The map util class for collision checking
 * @param Dim is the dimension of the workspace
 */
template <int Dim>
class MapUtil {
 public:
  // Constructor
  MapUtil(float res, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max,
          float inflation, float obst_max_vel) {
    /* --------- Initialize parameters --------- */
    setInflation(inflation);                               // Set inflation
    setResolution(res);                                    // Set the resolution
    setMapSize(x_min, x_max, y_min, y_max, z_min, z_max);  // Set the cells and z_boundaries
    setObstMaxVelocity(obst_max_vel);                      // Set obstacle maximum velocity
  }

  // Copy constructor (needed because std::mutex is not copyable)
  MapUtil(const MapUtil& other)
      : map_(other.map_),
        heat_(other.heat_),
        dynamic_heat_enabled_(other.dynamic_heat_enabled_),
        dynamic_as_occupied_current_(other.dynamic_as_occupied_current_),
        dynamic_as_occupied_future_(other.dynamic_as_occupied_future_),
        heat_w_(other.heat_w_),
        heat_alpha0_(other.heat_alpha0_),
        heat_alpha1_(other.heat_alpha1_),
        heat_p_(other.heat_p_),
        heat_q_(other.heat_q_),
        heat_tau_ratio_(other.heat_tau_ratio_),
        heat_gamma_(other.heat_gamma_),
        heat_Hmax_(other.heat_Hmax_),
        dyn_base_inflation_m_(other.dyn_base_inflation_m_),
        dyn_heat_tube_radius_m_(other.dyn_heat_tube_radius_m_),
        heat_num_samples_(other.heat_num_samples_),
        dyn_pred_samples_(other.dyn_pred_samples_),
        dyn_pred_times_(other.dyn_pred_times_),
        static_heat_enabled_(other.static_heat_enabled_),
        static_heat_alpha_(other.static_heat_alpha_),
        static_heat_p_(other.static_heat_p_),
        static_heat_Hmax_(other.static_heat_Hmax_),
        static_heat_rmax_m_(other.static_heat_rmax_m_),
        static_heat_default_radius_m_(other.static_heat_default_radius_m_),
        static_heat_boundary_only_(other.static_heat_boundary_only_),
        static_heat_apply_on_unknown_(other.static_heat_apply_on_unknown_),
        static_heat_exclude_dynamic_(other.static_heat_exclude_dynamic_),
        static_heat_radius_fn_(other.static_heat_radius_fn_),
        static_heat_off_(other.static_heat_off_),
        static_heat_off_Rcell_(other.static_heat_off_Rcell_),
        static_heat_off_res_(other.static_heat_off_res_),
        static_heat_off_rmax_m_(other.static_heat_off_rmax_m_),
        // static_heat_mutex_ is default-constructed (mutexes cannot be copied)
        use_soft_cost_obstacles_(other.use_soft_cost_obstacles_),
        obstacle_soft_cost_(other.obstacle_soft_cost_),
        heat_cutoff_ratio_(other.heat_cutoff_ratio_),
        // 2D ground robot map data
        map_2d_(other.map_2d_),
        heat_2d_(other.heat_2d_),
        terrain_cost_(other.terrain_cost_),
        col_min_height_(other.col_min_height_),
        col_max_height_(other.col_max_height_),
        has_2d_map_(other.has_2d_map_),
        res_(other.res_),
        total_size_(other.total_size_),
        inflation_(other.inflation_),
        origin_d_(other.origin_d_),
        center_map_(other.center_map_),
        dim_(other.dim_),
        prev_dim_(other.prev_dim_),
        dim_xy_(other.dim_xy_),
        x_map_min_(other.x_map_min_),
        x_map_max_(other.x_map_max_),
        y_map_min_(other.y_map_min_),
        y_map_max_(other.y_map_max_),
        z_map_min_(other.z_map_min_),
        z_map_max_(other.z_map_max_),
        x_min_(other.x_min_),
        x_max_(other.x_max_),
        y_min_(other.y_min_),
        y_max_(other.y_max_),
        z_min_(other.z_min_),
        z_max_(other.z_max_),
        obst_max_vel_(other.obst_max_vel_),
        cells_x_(other.cells_x_),
        cells_y_(other.cells_y_),
        cells_z_(other.cells_z_),
        val_occ_(other.val_occ_),
        val_free_(other.val_free_),
        val_unknown_(other.val_unknown_),
        map_initialized_(other.map_initialized_),
        min_point_(other.min_point_),
        max_point_(other.max_point_) {
    // Mutex is default-constructed
  }

  // Destructor
  ~MapUtil() {
    // Clear the map
    map_.clear();
  }

  /** @brief Build the voxel map from point clouds, inflating obstacles and computing heat maps.
   *  @param cloud Point cloud of occupied observations.
   *  @param unknown_cloud Point cloud of unknown-space observations (currently unused).
   *  @param cells_x Number of cells in x dimension.
   *  @param cells_y Number of cells in y dimension.
   *  @param cells_z Number of cells in z dimension.
   *  @param center_map Center of the map window in world coordinates.
   *  @param z_ground Minimum z boundary (ground level).
   *  @param z_max Maximum z boundary (ceiling).
   *  @param inflation Inflation radius in meters for obstacle dilation.
   *  @param obst_pos Current positions of dynamic obstacles.
   *  @param obst_bbox Bounding box half-extents of dynamic obstacles.
   *  @param traj_max_time Time horizon for dynamic obstacle inflation and heat computation.
   */
  void readMap(const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& cloud,
               const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& unknown_cloud, int cells_x,
               int cells_y, int cells_z, const Vec3f& center_map, double z_ground, double z_max,
               double inflation, const vec_Vecf<3>& obst_pos, const vec_Vecf<3>& obst_bbox,
               double traj_max_time) {
    (void)unknown_cloud;  // unknown-space soft costs removed
    // 1) Compute X/Y dims with inflation pad
    int pad = int(std::ceil(5.0 * inflation / res_));
    int dimX = cells_x + pad, dimY = cells_y + pad, dimZ = cells_z;

    // 2) Compute how many cells below/above center we keep,
    //    strictly within [z_ground, z_max]
    // Split dimZ into down/up halves. Use ceiling for `up` so dimZ==1 yields
    // down=0/up=1 (a single z-slice at center) instead of collapsing to 0.
    int down = dimZ / 2;
    int up = (dimZ + 1) / 2;
    // world coords of bottom slice:
    float bot = center_map.z() - down * res_;
    if (bot < z_ground) down = std::max(int(std::floor((center_map.z() - z_ground) / res_)), 0);
    // top slice:
    float top = center_map.z() + up * res_;
    if (top > z_max) up = std::max(int(std::floor((z_max - center_map.z()) / res_)), 1);
    dimZ = std::max(down + up, 1);

    // 3) Compute origin (global coords of cell (0,0,0)) and snap to resolution grid
    //    so that voxel boundaries align with the global mapper's grid
    Vec3f origin;
    origin.x() = std::floor((center_map.x() - (dimX * res_) / 2.0f) / res_) * res_;
    origin.y() = std::floor((center_map.y() - (dimY * res_) / 2.0f) / res_) * res_;
    origin.z() = center_map.z() - down * res_;
    // ensure origin.z >= z_ground and origin.z+dimZ*res <= z_max
    origin.z() = std::clamp(origin.z(), z_ground, z_max - dimZ * res_);

    // 4) Allocate map and fill as unknown
    size_t total = size_t(dimX) * dimY * dimZ;
    map_.assign(total, val_unknown_);

    // Optional: track which occupied voxels were introduced by dynamic obstacle hard-blocking,
    // so static heat can exclude them and avoid double counting.
    const bool need_dyn_mask = static_heat_enabled_ && static_heat_exclude_dynamic_;
    std::vector<uint8_t> dyn_occ_mask;
    if (need_dyn_mask) dyn_occ_mask.assign(total, 0);

    auto mark_dyn_occ = [&](size_t lin) {
      if (need_dyn_mask) dyn_occ_mask[lin] = 1;
    };

    // 5) Precompute inflation offsets
    int m = int(std::floor(inflation / res_));
    std::vector<Vec3i> offsets;
    offsets.reserve((2 * m + 1) * (2 * m + 1) * (2 * m + 1));
    for (int dx = -m; dx <= m; ++dx)
      for (int dy = -m; dy <= m; ++dy)
        for (int dz = -m; dz <= m; ++dz) offsets.emplace_back(dx, dy, dz);

    // 6) Helpers for indexing
    auto idx3 = [&](int x, int y, int z) {
      return size_t(x) + size_t(dimX) * y + size_t(dimX) * size_t(dimY) * z;
    };

// 7) Rasterize & inflate, skipping points outside grid bounds
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < cloud->points.size(); ++i) {
      const auto& P = cloud->points[i];
      if (P.z < z_ground || P.z > z_max) continue;

      // Compute grid indices and SKIP (not clamp) points outside the grid.
      // Clamping would project distant obstacles onto grid boundaries,
      // creating false occupied voxels that move with the sliding window.
      int xi = int(std::floor((P.x - origin.x()) / res_));
      int yi = int(std::floor((P.y - origin.y()) / res_));
      int zi = int(std::floor((P.z - origin.z()) / res_));
      if (xi < 0 || xi >= dimX || yi < 0 || yi >= dimY || zi < 0 || zi >= dimZ) continue;

      // mark occupied
      map_[idx3(xi, yi, zi)] = val_occ_;
      // inflate neighborhood
      for (auto& off : offsets) {
        int x2 = xi + off.x(), y2 = yi + off.y(), z2 = zi + off.z();
        if (x2 < 0 || x2 >= dimX || y2 < 0 || y2 >= dimY || z2 < 0 || z2 >= dimZ) continue;
        map_[idx3(x2, y2, z2)] = val_occ_;
      }
    }

    // 8) Dynamic obstacles:
    // Option A (legacy): mark reachable sphere as occupied (hard).
    // Option B (new): keep occupancy purely static/unknown and compute a soft heat map for global
    // planning.

    // 8a) Optionally hard-block the *current* dynamic obstacle footprint (hard constraint).
    //     This is complementary to heat-based planning:
    //       - current pose can be treated as occupied (planner must not go through it)
    //       - future motion is represented via heat_ (soft cost) unless dynamic_as_occupied_future_
    //       is true.
    //
    //     Now using bounding box dimensions instead of spherical radius.
    if (dynamic_as_occupied_current_) {
      const double base_inflation = std::max((double)dyn_base_inflation_m_, (double)inflation);
      if (base_inflation >= 0.0 && !obst_pos.empty()) {
        for (size_t k = 0; k < obst_pos.size(); ++k) {
          const auto& O = obst_pos[k];

          // Get bbox dimensions for this obstacle (default to small cube if not provided)
          double bbox_x = 0.8, bbox_y = 0.8, bbox_z = 0.8;
          if (k < obst_bbox.size()) {
            bbox_x = obst_bbox[k].x();
            bbox_y = obst_bbox[k].y();
            bbox_z = obst_bbox[k].z();
          }

          // Half-extents plus inflation
          const double hx = bbox_x + base_inflation;
          const double hy = bbox_y + base_inflation;
          const double hz = bbox_z + base_inflation;

          int ix_min = int(std::floor((O.x() - hx - origin.x()) / res_));
          int ix_max = int(std::floor((O.x() + hx - origin.x()) / res_));
          int iy_min = int(std::floor((O.y() - hy - origin.y()) / res_));
          int iy_max = int(std::floor((O.y() + hy - origin.y()) / res_));
          int iz_min = int(std::floor((O.z() - hz - origin.z()) / res_));
          int iz_max = int(std::floor((O.z() + hz - origin.z()) / res_));

          ix_min = std::clamp(ix_min, 0, dimX - 1);
          ix_max = std::clamp(ix_max, 0, dimX - 1);
          iy_min = std::clamp(iy_min, 0, dimY - 1);
          iy_max = std::clamp(iy_max, 0, dimY - 1);
          iz_min = std::clamp(iz_min, 0, dimZ - 1);
          iz_max = std::clamp(iz_max, 0, dimZ - 1);

          for (int ix = ix_min; ix <= ix_max; ++ix) {
            const double xc = origin.x() + (ix + 0.5) * res_;
            const double dx = std::abs(xc - O.x());
            for (int iy = iy_min; iy <= iy_max; ++iy) {
              const double yc = origin.y() + (iy + 0.5) * res_;
              const double dy = std::abs(yc - O.y());
              for (int iz = iz_min; iz <= iz_max; ++iz) {
                const double zc = origin.z() + (iz + 0.5) * res_;
                const double dz = std::abs(zc - O.z());

                // AABB-AABB intersection: check if cell volume intersects obstacle box
                // Cell occupies [center - res/2, center + res/2] in each dimension
                if (dx <= hx + res_ * 0.5 && dy <= hy + res_ * 0.5 && dz <= hz + res_ * 0.5) {
                  const size_t lin = idx3(ix, iy, iz);
                  map_[lin] = val_occ_;
                  mark_dyn_occ(lin);
                }
              }
            }
          }
        }
      }
    }  // end if (dynamic_as_occupied_current_)

    if (dynamic_as_occupied_future_) {
      const double motion_radius = (obst_max_vel_ * traj_max_time);  // [m] reachable distance
      if (motion_radius > 0.0 && !obst_pos.empty()) {
        for (size_t k = 0; k < obst_pos.size(); ++k) {
          const auto& O = obst_pos[k];

          // Get bbox dimensions for this obstacle (default to small cube if not provided)
          double bbox_x = 0.8, bbox_y = 0.8, bbox_z = 0.8;
          if (k < obst_bbox.size()) {
            bbox_x = obst_bbox[k].x();
            bbox_y = obst_bbox[k].y();
            bbox_z = obst_bbox[k].z();
          }

          // Reachable region: bbox half-extents plus motion radius
          const double hx = bbox_x + motion_radius;
          const double hy = bbox_y + motion_radius;
          const double hz = bbox_z + motion_radius;

          int ix_min = int(std::floor((O.x() - hx - origin.x()) / res_));
          int ix_max = int(std::floor((O.x() + hx - origin.x()) / res_));
          int iy_min = int(std::floor((O.y() - hy - origin.y()) / res_));
          int iy_max = int(std::floor((O.y() + hy - origin.y()) / res_));
          int iz_min = int(std::floor((O.z() - hz - origin.z()) / res_));
          int iz_max = int(std::floor((O.z() + hz - origin.z()) / res_));

          ix_min = std::clamp(ix_min, 0, dimX - 1);
          ix_max = std::clamp(ix_max, 0, dimX - 1);
          iy_min = std::clamp(iy_min, 0, dimY - 1);
          iy_max = std::clamp(iy_max, 0, dimY - 1);
          iz_min = std::clamp(iz_min, 0, dimZ - 1);
          iz_max = std::clamp(iz_max, 0, dimZ - 1);

          for (int ix = ix_min; ix <= ix_max; ++ix) {
            const double xc = origin.x() + (ix + 0.5) * res_;
            const double dx = std::abs(xc - O.x());

            for (int iy = iy_min; iy <= iy_max; ++iy) {
              const double yc = origin.y() + (iy + 0.5) * res_;
              const double dy = std::abs(yc - O.y());

              for (int iz = iz_min; iz <= iz_max; ++iz) {
                const double zc = origin.z() + (iz + 0.5) * res_;
                const double dz = std::abs(zc - O.z());

                // AABB-AABB intersection: check if cell volume intersects reachable box
                // Cell occupies [center - res/2, center + res/2] in each dimension
                if (dx <= hx + res_ * 0.5 && dy <= hy + res_ * 0.5 && dz <= hz + res_ * 0.5) {
                  const size_t lin = idx3(ix, iy, iz);
                  map_[lin] = val_occ_;
                  mark_dyn_occ(lin);
                }
              }
            }
          }
        }
      }
    }

    // 8b) Mark boundary walls — keep planner away from global y limits
    //     so SFC corridors have room to expand.
    //     Note: z boundaries are NOT marked because z_min/z_max already
    //     constrain the map. Marking them creates false heat sources that
    //     push the global planner away from the floor/ceiling.
    {
      const float buf = inflation;  // reuse the inflation radius as boundary buffer
      for (int ix = 0; ix < dimX; ++ix) {
        for (int iy = 0; iy < dimY; ++iy) {
          const float wy = origin.y() + (iy + 0.5f) * res_;
          if (wy <= y_map_min_ + buf || wy >= y_map_max_ - buf) {
            for (int iz = 0; iz < dimZ; ++iz) {
              map_[idx3(ix, iy, iz)] = val_occ_;
            }
          }
        }
      }
    }

    // 9) Update metadata
    dim_ = Veci<3>(dimX, dimY, dimZ);
    dim_xy_ = dimX * dimY;  // Cache the stride for fast indexing
    total_size_ = total;
    origin_d_ = origin;
    center_map_ = center_map;
    // 10) Build soft cost map (dynamic heat + unknown-space cost).
    //     This stays separate from hard occupancy (map_) and is consumed by weighted A*:
    //       edge_cost += heat_w_ * getHeat(...)
    const bool need_heat = dynamic_heat_enabled_ || static_heat_enabled_;
    if (need_heat) {
      heat_.assign(total_size_, 0.0f);
    } else {
      heat_.clear();
    }

    // 10a) Dynamic obstacle heat (time-invariant, max-over-time tube).
    //      Heat is used only by weighted A* (global planner); static occupancy remains hard.
    if (dynamic_heat_enabled_) {
      const float Th = std::max(0.0f, (float)traj_max_time);
      const float tau_w = std::max(1e-3f, heat_tau_ratio_ * std::max(1e-3f, Th));

      // Determine time samples:
      std::vector<float> t_samples;
      if (!dyn_pred_times_.empty()) {
        t_samples = dyn_pred_times_;
      } else {
        const int M = std::max(2, heat_num_samples_);
        t_samples.resize(M);
        for (int j = 0; j < M; ++j) t_samples[j] = (float)j * Th / (float)(M - 1);
      }

      const size_t K = obst_pos.size();
      const size_t J_pred = dyn_pred_samples_.size();
      const size_t T_pred = dyn_pred_times_.size();
      if (K > 0) {
        // Lightweight pow for small integer exponents.
        auto pow_fast = [](float x, int p) -> float {
          x = std::max(0.0f, x);
          switch (p) {
            case 1:
              return x;
            case 2:
              return x * x;
            case 3:
              return x * x * x;
            case 4: {
              const float x2 = x * x;
              return x2 * x2;
            }
            default:
              return std::pow(x, (float)p);
          }
        };

        // Precompute obstacle centers (float), bbox half-extents, and reachable radii.
        std::vector<Eigen::Vector3f> ck_list(K);
        std::vector<Eigen::Vector3f> hk_list(K);  // bbox half-extents
        std::vector<float> Rreach_list(K);

        for (size_t k = 0; k < K; ++k) {
          ck_list[k] = obst_pos[k].cast<float>();

          // Get bbox half-extents for this obstacle (default to small cube if not provided)
          float hx = 0.4f, hy = 0.4f, hz = 0.4f;  // default half-extents
          if (k < obst_bbox.size()) {
            hx = obst_bbox[k].x();
            hy = obst_bbox[k].y();
            hz = obst_bbox[k].z();
          }
          hk_list[k] = Eigen::Vector3f(hx, hy, hz);

          // Reachable radius: bbox extent + motion
          const float max_extent = std::max({hx, hy, hz});
          Rreach_list[k] = max_extent + (float)obst_max_vel_ * Th;
        }

        // Tube radius and time-decay weight per sample
        const float R0 = std::max(0.0f, dyn_heat_tube_radius_m_);
        const size_t J = t_samples.size();
        std::vector<float> Rj(J), Wj(J);
        for (size_t j = 0; j < J; ++j) {
          const float tj = std::max(0.0f, t_samples[j]);
          Rj[j] = R0 + heat_gamma_ * tj;
          Wj[j] = std::exp(-tj / tau_w);
        }

        // Precompute predicted centers (fallback to ck if unavailable).
        std::vector<Eigen::Vector3f> cj_flat(K * J);
        for (size_t k = 0; k < K; ++k) {
          for (size_t j = 0; j < J; ++j) {
            Eigen::Vector3f cj = ck_list[k];
            if (k < dyn_pred_samples_.size() && j < dyn_pred_samples_[k].size())
              cj = dyn_pred_samples_[k][j].cast<float>();
            cj_flat[k * J + j] = cj;
          }
        }

        const int dim0 = dim_(0);
        const int dim1 = dim_(1);
        const int plane = dim0 * dim1;

#pragma omp parallel for schedule(static)
        for (int idx = 0; idx < total_size_; ++idx) {
          // Heat is only relevant for traversable cells; skip hard obstacles unless soft-cost mode
          if (map_[idx] > val_free_ && !use_soft_cost_obstacles_) continue;

          const int ix = idx % dim0;
          const int iy = (idx / dim0) % dim1;
          const int iz = idx / plane;

          const float xw = origin_d_.x() + (ix + 0.5f) * res_;
          const float yw = origin_d_.y() + (iy + 0.5f) * res_;
          const float zw = origin_d_.z() + (iz + 0.5f) * res_;

          float best = 0.0f;

          for (size_t k = 0; k < K; ++k) {
            // Base reachable radius (finite horizon)
            float Hbase = 0.0f;
            const float Rreach = Rreach_list[k];
            if (Rreach > 1e-6f) {
              const Eigen::Vector3f& ck = ck_list[k];
              const Eigen::Vector3f& hk = hk_list[k];

              // Compute distance from point to box (0 if inside, positive if outside)
              const float dx_abs = std::abs(xw - ck.x());
              const float dy_abs = std::abs(yw - ck.y());
              const float dz_abs = std::abs(zw - ck.z());

              const float dx_box = std::max(0.0f, dx_abs - hk.x());
              const float dy_box = std::max(0.0f, dy_abs - hk.y());
              const float dz_box = std::max(0.0f, dz_abs - hk.z());

              const float d2 = dx_box * dx_box + dy_box * dy_box + dz_box * dz_box;
              const float R2 = Rreach * Rreach;

              if (d2 <= R2) {
                const float d = std::sqrt(std::max(0.0f, d2));
                const float u = std::min(1.0f, std::max(0.0f, d / Rreach));
                Hbase = heat_alpha0_ * pow_fast(1.0f - u, heat_p_);
              }
            }

            // Tube bonus (max over time), bbox + growing margin with time, time-decayed
            float tube_max = 0.0f;
            const Eigen::Vector3f* cj_ptr = &cj_flat[k * J];
            const Eigen::Vector3f& hk = hk_list[k];

            for (size_t j = 0; j < J; ++j) {
              const float R = Rj[j];
              if (R <= 1e-6f) continue;

              const Eigen::Vector3f& cj = cj_ptr[j];

              // Distance from point to box at predicted position
              const float dx_abs = std::abs(xw - cj.x());
              const float dy_abs = std::abs(yw - cj.y());
              const float dz_abs = std::abs(zw - cj.z());

              const float dx_box = std::max(0.0f, dx_abs - hk.x());
              const float dy_box = std::max(0.0f, dy_abs - hk.y());
              const float dz_box = std::max(0.0f, dz_abs - hk.z());

              const float d2 = dx_box * dx_box + dy_box * dy_box + dz_box * dz_box;
              const float R2 = R * R;

              if (d2 > R2) continue;

              const float d = std::sqrt(std::max(0.0f, d2));
              const float u = std::min(1.0f, std::max(0.0f, d / R));
              const float g = pow_fast(1.0f - u, heat_q_);
              tube_max = std::max(tube_max, Wj[j] * g);
            }

            float Hk = Hbase + heat_alpha1_ * tube_max;
            if (heat_Hmax_ > 0.0f) Hk = std::min(Hk, heat_Hmax_);

            best = std::max(best, Hk);
          }

          heat_[idx] = std::max(heat_[idx], best);
        }
      }
    }

    // 10c) Static obstacle heat (soft cost halo around occupied voxels).
    // Performance strategy:
    //  - Use only boundary occupied voxels as seeds (6-neighborhood boundary test).
    //  - Apply a radial falloff halo in FREE (and optionally UNKNOWN) space.
    //  - Radius is provided by a user function (lambda) evaluated at the seed voxel center.
    if (static_heat_enabled_ && static_heat_alpha_ > 0.0f && static_heat_rmax_m_ > 1e-6f) {
      // Precompute offsets up to rmax once (distance stored in meters).
      const int Rcell = int(std::ceil(static_heat_rmax_m_ / res_));

      ensureStaticHeatOffsets(Rcell);
      const auto& off = static_heat_off_;

      auto idx3_local = [&](int x, int y, int z) {
        return size_t(x) + size_t(dimX) * size_t(y) + size_t(dimX) * size_t(dimY) * size_t(z);
      };

      // Collect seeds (boundary occupied voxels).
      std::vector<size_t> seeds;
      seeds.reserve(size_t(total_size_ / 50) + 1);

      const int nx[6] = {+1, -1, 0, 0, 0, 0};
      const int ny[6] = {0, 0, +1, -1, 0, 0};
      const int nz[6] = {0, 0, 0, 0, +1, -1};

      for (int z = 0; z < dimZ; ++z) {
        for (int y = 0; y < dimY; ++y) {
          for (int x = 0; x < dimX; ++x) {
            const size_t lin = idx3_local(x, y, z);
            if (map_[lin] != val_occ_) continue;
            if (need_dyn_mask && dyn_occ_mask[lin]) continue;

            if (!static_heat_boundary_only_) {
              seeds.push_back(lin);
              continue;
            }

            bool boundary = false;
            for (int k = 0; k < 6; ++k) {
              const int x2 = x + nx[k], y2 = y + ny[k], z2 = z + nz[k];
              if (x2 < 0 || x2 >= dimX || y2 < 0 || y2 >= dimY || z2 < 0 || z2 >= dimZ) {
                boundary = true;
                break;
              }
              const size_t nlin = idx3_local(x2, y2, z2);
              if (map_[nlin] != val_occ_) {
                boundary = true;
                break;
              }
            }

            if (boundary) seeds.push_back(lin);
          }
        }
      }

      // Apply halo (max aggregation) around each seed
      for (const auto& lin : seeds) {
        const int x0 = int(lin % size_t(dimX));
        const int y0 = int((lin / size_t(dimX)) % size_t(dimY));
        const int z0 = int(lin / (size_t(dimX) * size_t(dimY)));

        // Seed voxel center in world coordinates (for radius function)
        const float xc0 = origin_d_.x() + (x0 + 0.5f) * res_;
        const float yc0 = origin_d_.y() + (y0 + 0.5f) * res_;
        const float zc0 = origin_d_.z() + (z0 + 0.5f) * res_;
        const Eigen::Vector3f seed_world(xc0, yc0, zc0);

        float Rm = static_heat_default_radius_m_;
        if (static_heat_radius_fn_) Rm = static_heat_radius_fn_(seed_world);

        // Clamp radius for safety/perf
        Rm = std::clamp(Rm, 0.0f, static_heat_rmax_m_);
        if (Rm <= 1e-6f) continue;

        for (const auto& o : off) {
          if (o.d_m > Rm) continue;

          const int x = x0 + o.dx, y = y0 + o.dy, z = z0 + o.dz;
          if (x < 0 || x >= dimX || y < 0 || y >= dimY || z < 0 || z >= dimZ) continue;

          const size_t idx = idx3_local(x, y, z);

          const bool target_occupied = map_[idx] > val_free_;

          // Skip occupied cells unless soft-cost mode is on
          if (target_occupied && !use_soft_cost_obstacles_) continue;

          float w;
          if (target_occupied) {
            // Occupied cells get max heat (they ARE the obstacle)
            w = static_heat_alpha_;
          } else {
            // Free cells get distance-decay
            const float u = std::min(1.0f, std::max(0.0f, o.d_m / Rm));
            const float base = 1.0f - u;
            float power_result;
            if (static_heat_p_ == 2) {
              power_result = base * base;
            } else if (static_heat_p_ == 3) {
              power_result = base * base * base;
            } else if (static_heat_p_ == 4) {
              const float base2 = base * base;
              power_result = base2 * base2;
            } else {
              power_result = std::pow(base, float(static_heat_p_));
            }
            w = static_heat_alpha_ * power_result;
          }

          if (static_heat_Hmax_ > 0.0f) w = std::min(w, static_heat_Hmax_);

          if (w > 0.0f) {
            heat_[idx] = std::max(heat_[idx], w);
          }
        }
      }
    }
  }

  // ---------------- Dynamic heat-map API ----------------

  /** @brief Enable or disable hard-blocking of the current dynamic obstacle position. */
  void setDynamicAsOccupiedCurrentPos(bool enabled) { dynamic_as_occupied_current_ = enabled; }

  /** @brief Enable or disable hard-blocking of future reachable dynamic obstacle positions. */
  void setDynamicAsOccupiedFuturePos(bool enabled) { dynamic_as_occupied_future_ = enabled; }

  /** @brief Enable or disable dynamic obstacle heat map computation. */
  void setDynamicHeatEnabled(bool enabled) { dynamic_heat_enabled_ = enabled; }

  /** @brief Set the global planner weight for heat cost (edge_cost += w_heat * heat). */
  void setHeatWeight(float w_heat) { heat_w_ = w_heat; }

  /** @brief Get the current heat weight used by the global planner.
   *  @return Heat weight value.
   */
  float getHeatWeight() const { return heat_w_; }

  /** @brief Check whether dynamic heat map computation is enabled.
   *  @return True if dynamic heat is enabled.
   */
  bool dynamicHeatEnabled() const { return dynamic_heat_enabled_; }

  // ---------------- Static obstacle heat-map API ----------------

  /** @brief Enable or disable static obstacle heat halo accumulation. */
  void setStaticHeatEnabled(bool enabled) { static_heat_enabled_ = enabled; }

  /** @brief Check whether static obstacle heat computation is enabled.
   *  @return True if static heat is enabled.
   */
  bool staticHeatEnabled() const { return static_heat_enabled_; }

  /** @brief Configure static heat map shaping parameters.
   *  @param alpha Peak heat magnitude at distance zero.
   *  @param p Falloff power exponent in (1 - d/R)^p.
   *  @param Hmax Maximum heat cap (<=0 disables capping).
   *  @param rmax_m Maximum halo radius in meters.
   *  @param boundary_only If true, use only boundary occupied voxels as seeds.
   *  @param apply_on_unknown If true, apply halo to unknown cells in addition to free cells.
   *  @param exclude_dynamic If true, exclude dynamically occupied voxels from static heat seeds.
   */
  void setStaticHeatParams(float alpha, int p, float Hmax, float rmax_m, bool boundary_only = true,
                           bool apply_on_unknown = false, bool exclude_dynamic = true) {
    static_heat_alpha_ = std::max(0.0f, alpha);
    static_heat_p_ = std::max(1, p);
    static_heat_Hmax_ = Hmax;
    static_heat_rmax_m_ = std::max(0.0f, rmax_m);
    static_heat_boundary_only_ = boundary_only;
    static_heat_apply_on_unknown_ = apply_on_unknown;
    static_heat_exclude_dynamic_ = exclude_dynamic;
  }

  /** @brief Set a per-seed radius function for static heat halo computation.
   *  @param fn Function mapping seed voxel center (world coords) to halo radius; empty uses
   * default.
   *  @param default_radius_m Fallback radius when fn is not set.
   */
  void setStaticHeatRadiusFunction(const std::function<float(const Eigen::Vector3f&)>& fn,
                                   float default_radius_m) {
    static_heat_radius_fn_ = fn;
    static_heat_default_radius_m_ = std::max(0.0f, default_radius_m);
  }

  /** @brief Configure dynamic heat map shaping parameters.
   *  @param alpha0 Base heat magnitude from reachable sphere.
   *  @param alpha1 Tube bonus heat magnitude from predicted trajectory.
   *  @param p Falloff power for reachable sphere heat.
   *  @param q Falloff power for tube bonus heat.
   *  @param tau_w_ratio Time decay ratio (tau_w = tau_w_ratio * T_horizon).
   *  @param gamma Tube radius growth rate in meters per second.
   *  @param Hmax Maximum heat cap.
   *  @param base_inflation_m Base inflation radius for current obstacle position.
   */
  void setDynamicHeatParams(float alpha0, float alpha1, int p, int q, float tau_w_ratio,
                            float gamma, float Hmax, float base_inflation_m) {
    heat_alpha0_ = alpha0;
    heat_alpha1_ = alpha1;
    heat_p_ = std::max(1, p);
    heat_q_ = std::max(1, q);
    heat_tau_ratio_ = std::max(1e-3f, tau_w_ratio);
    heat_gamma_ = std::max(0.0f, gamma);
    heat_Hmax_ = std::max(0.0f, Hmax);
    dyn_base_inflation_m_ = std::max(0.0f, base_inflation_m);
  }

  /** @brief Set the base tube radius for dynamic heat corridor computation. */
  void setDynHeatTubeRadius(float r) { dyn_heat_tube_radius_m_ = r; }

  /** @brief Enable soft-cost mode where occupied cells receive a finite cost instead of being
   * blocked.
   *  @param enable If true, occupied cells are traversable with a soft cost penalty.
   *  @param cost The soft cost value assigned to occupied cells.
   */
  void setSoftCostObstacles(bool enable, float cost) {
    use_soft_cost_obstacles_ = enable;
    obstacle_soft_cost_ = cost;
  }
  /** @brief Check whether soft-cost obstacle mode is enabled.
   *  @return True if soft-cost mode is active.
   */
  bool useSoftCostObstacles() const { return use_soft_cost_obstacles_; }

  /** @brief Get the soft cost value assigned to occupied cells.
   *  @return Soft cost value.
   */
  float getObstacleSoftCost() const { return obstacle_soft_cost_; }

  /** @brief Set predicted trajectory samples for dynamic obstacle heat tube computation.
   *  @param pred_samples Predicted positions per obstacle, indexed [obstacle][time_step].
   *  @param pred_times Time stamps aligned with the prediction samples.
   */
  void setDynamicPredictedSamples(const std::vector<vec_Vecf<3>>& pred_samples,
                                  const std::vector<float>& pred_times) {
    dyn_pred_samples_ = pred_samples;
    dyn_pred_times_ = pred_times;
  }

  /** @brief Query the heat value at integer voxel coordinates.
   *  @param x Voxel x index.
   *  @param y Voxel y index.
   *  @param z Voxel z index.
   *  @return Heat value at the given voxel, or 0 if out of bounds or heat map is empty.
   */
  float getHeat(int x, int y, int z) const {
    if (heat_.empty()) return 0.0f;
    if (x < 0 || x >= dim_(0) || y < 0 || y >= dim_(1) || z < 0 || z >= dim_(2)) return 0.0f;
    const int idx = x + y * dim_(0) + z * dim_(0) * dim_(1);
    if (idx < 0 || idx >= (int)heat_.size()) return 0.0f;
    return heat_[(size_t)idx];
  }

  /** @brief Get 2D heat value at grid (x,y). Falls back to 3D heat at z=0 if no 2D heat. */
  float getHeat2D(int x, int y) const {
    if (!heat_2d_.empty()) {
      if (x < 0 || x >= dim_(0) || y < 0 || y >= dim_(1)) return 0.0f;
      const size_t idx = static_cast<size_t>(x) + static_cast<size_t>(dim_(0)) * y;
      if (idx >= heat_2d_.size()) return 0.0f;
      return heat_2d_[idx];
    }
    return getHeat(x, y, 0);
  }

  /** @brief Get positions of voxels with heat above a threshold for visualization.
   *  @param threshold Minimum heat value to include.
   *  @return Vector of world-coordinate positions of heat voxels.
   */
  vec_Vecf<3> getHeatCloud(float threshold = 0.01f) const {
    vec_Vecf<3> cloud;
    if (heat_.empty()) return cloud;

    for (int x = 0; x < dim_(0); ++x) {
      for (int y = 0; y < dim_(1); ++y) {
        for (int z = 0; z < dim_(2); ++z) {
          const int idx = x + y * dim_(0) + z * dim_xy_;
          if (idx >= 0 && idx < (int)heat_.size() && heat_[idx] > threshold) {
            Vecf<3> pt;
            pt(0) = origin_d_(0) + (x + 0.5f) * res_;
            pt(1) = origin_d_(1) + (y + 0.5f) * res_;
            pt(2) = origin_d_(2) + (z + 0.5f) * res_;
            cloud.push_back(pt);
          }
        }
      }
    }
    return cloud;
  }

  /** @brief Get the raw heat values for all voxels.
   *  @return Copy of the heat value vector.
   */
  std::vector<float> getHeatValues() const { return heat_; }

  /** @brief Get the maximum heat value across all voxels.
   *  @return Maximum heat value, or 0 if the heat map is empty.
   */
  float getMaxHeat() const {
    if (heat_.empty()) return 0.0f;
    return *std::max_element(heat_.begin(), heat_.end());
  }

  /** @brief Precompute 3D grid offsets for a given inflation radius in cells.
   *  @param inflation_cells Number of cells to inflate in each dimension.
   *  @return Vector of integer offset vectors within the inflation cube.
   */
  vec_Veci<3> computeInflationOffsets(const Veci<3>& inflation_cells) {
    vec_Veci<3> offsets;

    // include diagonal offsets
    for (int dx = -inflation_cells(0); dx <= inflation_cells(0); ++dx) {
      for (int dy = -inflation_cells(1); dy <= inflation_cells(1); ++dy) {
        for (int dz = -inflation_cells(2); dz <= inflation_cells(2); ++dz) {
          offsets.push_back(Veci<3>(dx, dy, dz));
        }
      }
    }

    return offsets;
  }

  /** @brief Convert a linear map index to 3D voxel coordinates.
   *  @param index Linear index into the map array.
   *  @return 3D integer voxel coordinates.
   */
  Veci<3> indexToVeci3(int index) {
    Veci<3> position;
    position[0] = index % dim_(0);
    position[1] = (index / dim_(0)) % dim_(1);
    position[2] = index / (dim_(0) * dim_(1));
    return position;
  }

  /** @brief Set the number of cells in each dimension. */
  void setCellSize(int cells_x, int cells_y, int cells_z) {
    // Set cells
    cells_x_ = cells_x;
    cells_y_ = cells_y;
    cells_z_ = cells_z;
  }

  /** @brief Set the world-coordinate boundaries of the map. */
  void setMapSize(float x_min, float x_max, float y_min, float y_max, float z_min, float z_max) {
    // Set map boundaries
    x_map_min_ = x_min;
    x_map_max_ = x_max;
    y_map_min_ = y_min;
    y_map_max_ = y_max;
    z_map_min_ = z_min;
    z_map_max_ = z_max;
  }

  /** @brief Set the assumed maximum velocity of dynamic obstacles for inflation radius computation.
   */
  void setObstMaxVelocity(float obst_max_vel) {
    // Set obstacle maximum velocity
    obst_max_vel_ = obst_max_vel;
  }

  /**
   * @brief  Find a free point in the map that is closest to the given point
   * @param  vec_Vecf<3> point : The given point
   * @param  vec_Vecf<3> free_point : The free point that is closest to the given point
   * @return void
   */
  void findClosestFreePoint(const Vec3f& point, Vec3f& closest_free_point) {
    // Initialize the closest free point
    closest_free_point = point;

    // Check if the map is initialized
    if (!map_initialized_) {
      std::cout << "Map is not initialized" << std::endl;
      return;
    }

    // Get the position of the point in int
    Veci<3> point_int = floatToInt(point);

    // Get the index
    int index = getIndex(point_int);

    if (index >= 0 && index < total_size_) {
      // Check if the point is free
      if (map_[index] == val_free_) {
        closest_free_point = point;
        return;
      }

      // Get the neighboring indices
      std::vector<int> neighbor_indices;

      // Increase the radius until a free point is found
      for (float radius = 1.0; radius < 5.0;
           radius += 0.5)
      {
        neighbor_indices.clear();
        getNeighborIndices(point_int, neighbor_indices, radius);

        // Find the closest free point
        float min_dist = std::numeric_limits<float>::max();
        for (int neighbor_index : neighbor_indices) {
          if (neighbor_index >= 0 && neighbor_index < total_size_) {
            if (map_[neighbor_index] == val_free_) {
              Veci<3> neighbor_int = indexToVeci3(neighbor_index);
              Vec3f neighbor = intToFloat(neighbor_int);
              float dist = (neighbor - point).norm();
              if (dist < min_dist) {
                min_dist = dist;
                closest_free_point = neighbor;
              }
            }
          }
        }

        // Check if a free point is found
        if (min_dist < std::numeric_limits<float>::max()) {
          return;
        }
      }
    }
  }

  /**
   * @brief  Find a non-occupied point (free or unknown) in the map closest to
   *         the given point. Mirrors findClosestFreePoint but the predicate is
   *         "not occupied" instead of "free", so unknown cells are accepted.
   * @param  point : The given query position in world coordinates
   * @param  closest_non_occupied_point : Output position of the nearest non-occupied voxel
   */
  void findClosestNonOccupiedPoint(const Vec3f& point, Vec3f& closest_non_occupied_point) {
    closest_non_occupied_point = point;

    if (!map_initialized_) {
      std::cout << "Map is not initialized" << std::endl;
      return;
    }

    Veci<3> point_int = floatToInt(point);
    int index = getIndex(point_int);

    if (index >= 0 && index < total_size_) {
      // Already non-occupied (free or unknown): nothing to do.
      if (map_[index] <= val_free_) {
        closest_non_occupied_point = point;
        return;
      }

      std::vector<int> neighbor_indices;

      // Expanding-radius search until a non-occupied cell is found.
      for (float radius = 1.0; radius < 5.0; radius += 0.5) {
        neighbor_indices.clear();
        getNeighborIndices(point_int, neighbor_indices, radius);

        float min_dist = std::numeric_limits<float>::max();
        for (int neighbor_index : neighbor_indices) {
          if (neighbor_index >= 0 && neighbor_index < total_size_) {
            // val_unknown_ = -1, val_free_ = 0, occupied > 0 → "<= val_free_" is non-occupied.
            if (map_[neighbor_index] <= val_free_) {
              Veci<3> neighbor_int = indexToVeci3(neighbor_index);
              Vec3f neighbor = intToFloat(neighbor_int);
              float dist = (neighbor - point).norm();
              if (dist < min_dist) {
                min_dist = dist;
                closest_non_occupied_point = neighbor;
              }
            }
          }
        }

        if (min_dist < std::numeric_limits<float>::max()) {
          return;
        }
      }
    }
  }

  /**
   * @brief Get indices of the neighbors of a point given the radius
   * @param Veci<3> point_int : The given point
   * @param std::vector<int>& neighbor_indices : The indices of the neighbors
   * @param float radius : The radius
   * @return void
   * */
  void getNeighborIndices(const Veci<3>& point_int, std::vector<int>& neighbor_indices,
                          float radius) {
    // Get the radius in int
    float radius_int = radius / res_;
    Veci<3> radius_int_vec(radius_int, radius_int, radius_int);

    // Get the min and max positions
    Veci<3> min_pos = point_int - radius_int_vec;
    Veci<3> max_pos = point_int + radius_int_vec;

    // Iterate over the neighbors
    for (int x = min_pos[0]; x <= max_pos[0]; ++x) {
      for (int y = min_pos[1]; y <= max_pos[1]; ++y) {
        for (int z = min_pos[2]; z <= max_pos[2]; ++z) {
          // Check if the neighbor is inside the map
          if (x >= 0 && x < dim_(0) && y >= 0 && y < dim_(1) && z >= 0 && z < dim_(2)) {
            Veci<3> neighbor_int(x, y, z);
            int index = getIndex(neighbor_int);
            if (index >= 0 && index < total_size_) {
              neighbor_indices.push_back(index);
            }
          }
        }
      }
    }
  }

  // Check if the given point has any occupied neighbors.
  // Returns true if at least one neighboring cell is non-free.
  inline bool checkIfPointHasNonFreeNeighbour(const Veci<Dim>& pt) const {
    if constexpr (Dim == 2) {
      for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          // Skip the center point
          if (dx == 0 && dy == 0) continue;
          Veci<2> neighbor = pt;
          neighbor(0) += dx;
          neighbor(1) += dy;
          // Check if the neighbor is within the map and non-free.
          if (!isOutside(neighbor) && !isFree(neighbor)) return true;
        }
      }
    } else if constexpr (Dim == 3) {
      for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dz = -1; dz <= 1; ++dz) {
            // Skip the center point
            if (dx == 0 && dy == 0 && dz == 0) continue;
            Veci<3> neighbor = pt;
            neighbor(0) += dx;
            neighbor(1) += dy;
            neighbor(2) += dz;
            // Check if the neighbor is within the map and non-free.
            if (!isOutside(neighbor) && !isFree(neighbor)) return true;
          }
        }
      }
    }
    return false;
  }

  /** @brief Set the obstacle inflation radius in meters. */
  void setInflation(float inflation) { inflation_ = inflation; }

  /** @brief Set the voxel grid resolution in meters. */
  void setResolution(float res) { res_ = res; }

  /** @brief Check if a point lies inside an axis-aligned box defined by 8 corner vertices.
   *  @param p Query point.
   *  @param vertices Eight corner points of the axis-aligned box.
   *  @return True if the point is within the box (with a 0.5m margin).
   */
  bool isPointInBox(const Vec3f& p, const std::vector<Vec3f>& vertices) {
    // Assuming vertices contains 8 points representing the corners of the box in 3D space.
    // The box is axis-aligned and the vertices are provided in any order.

    // Find the minimum and maximum x, y, and z coordinates among the vertices.
    Vec3f min = vertices[0];
    Vec3f max = vertices[0];

    for (const auto& vertex : vertices) {
      min.x() = std::min(min.x(), vertex.x());
      min.y() = std::min(min.y(), vertex.y());
      min.z() = std::min(min.z(), vertex.z());

      max.x() = std::max(max.x(), vertex.x());
      max.y() = std::max(max.y(), vertex.y());
      max.z() = std::max(max.z(), vertex.z());
    }

    // Add buffer to the min and max coordinates to give some margin around the box.
    min -= Vec3f(0.5, 0.5, 0.5);  // min -= Vec3f
    max += Vec3f(0.5, 0.5, 0.5);  // max += Vec3f

    // Check if the point lies within the bounds defined by the min and max coordinates.
    return (p.x() >= min.x() && p.x() <= max.x()) && (p.y() >= min.y() && p.y() <= max.y()) &&
           (p.z() >= min.z() && p.z() <= max.z());
  }

  /** @brief Get the voxel grid resolution in meters.
   *  @return Resolution value.
   */
  decimal_t getRes() { return res_; }

  /** @brief Get the grid dimensions in cells.
   *  @return Integer vector of grid dimensions.
   */
  Veci<Dim> getDim() { return dim_; }

  /** @brief Get the world-coordinate origin of the grid (cell 0,0,0 corner).
   *  @return Float vector of origin coordinates.
   */
  Vecf<Dim> getOrigin() { return origin_d_; }

  /** @brief Convert voxel integer coordinates to a linear array index.
   *  @param pn Integer voxel coordinates.
   *  @return Linear index into the map array.
   */
  inline int getIndex(const Veci<Dim>& pn) const {
    return Dim == 2 ? pn(0) + dim_(0) * pn(1) : pn(0) + dim_(0) * pn(1) + dim_(0) * dim_(1) * pn(2);
  }
  // Get index of a cell in old map
  inline int getOldIndex(const Veci<Dim>& pn) const {
    return Dim == 2 ? pn(0) + prev_dim_(0) * pn(1)
                    : pn(0) + prev_dim_(0) * pn(1) + prev_dim_(0) * prev_dim_(1) * pn(2);
  }

  /** @brief Convert a linear map index to Dim-dimensional voxel coordinates.
   *  @param idx Linear index into the map array.
   *  @return Integer voxel coordinates.
   */
  Veci<Dim> getVoxelPos(int idx) {
    Veci<Dim> pn;
    if (Dim == 2) {
      pn(0) = idx % dim_(0);
      pn(1) = idx / dim_(0);
    } else {
      pn(0) = idx % dim_(0);
      pn(1) = (idx / dim_(0)) % dim_(1);
      pn(2) = idx / (dim_(0) * dim_(1));
    }
    return pn;
  }
  // Check if the given cell is outside of the map in i-the dimension
  inline bool isOutsideXYZ(const Veci<Dim>& n, int i) const { return n(i) < 0 || n(i) >= dim_(i); }
  // Check if the cell is free by index
  inline bool isFree(int idx) const { return map_[idx] == val_free_; }
  // Check if the cell is unknown by index
  inline bool isUnknown(int idx) const { return map_[idx] == val_unknown_; }
  // Check if the cell is occupied by index
  inline bool isOccupied(int idx) const { return map_[idx] > val_free_; }

  /** @brief Mark a voxel as occupied by integer coordinates. */
  inline void setOccupied(const Veci<Dim>& pn) {
    int index = getIndex(pn);
    if (index >= 0 && index < total_size_) {  // check that the point is inside the map
      map_[getIndex(pn)] = val_occ_;
    }
  }

  /** @brief Mark a voxel as free by integer coordinates. */
  inline void setFree(const Veci<Dim>& pn) {
    int index = getIndex(pn);
    if (index >= 0 && index < total_size_) {  // check that the point is inside the map
      map_[index] = val_free_;
    }
  }

  /** @brief Set all voxels within a cube of half-side d around center to free.
   *  @param center Integer voxel coordinates of the cube center.
   *  @param d Half-side length in meters (converted to voxel count internally).
   */
  inline void setFreeVoxelAndSurroundings(const Veci<Dim>& center, const float d) {
    int n_voxels = std::round(d / res_ + 0.5);  // convert distance to number of voxels
    for (int ix = -n_voxels; ix <= n_voxels; ix++) {
      for (int iy = -n_voxels; iy <= n_voxels; iy++) {
        for (int iz = -n_voxels; iz <= n_voxels; iz++) {
          Veci<Dim> voxel =
              center + Veci<Dim>(ix, iy, iz);  // Int coordinates of the voxel I'm going to clear

          setFree(voxel);
        }
      }
    }
  }

  // Check if the cell is outside by coordinate
  inline bool isOutsideOldMap(const Veci<Dim>& pn) const {
    for (int i = 0; i < Dim; i++)
      if (pn(i) < 0 || pn(i) >= prev_dim_(i)) return true;
    return false;
  }
  // Check if the cell is outside by coordinate
  inline bool isOutside(const Veci<Dim>& pn) const {
    for (int i = 0; i < Dim; i++)
      if (pn(i) < 0 || pn(i) >= dim_(i)) return true;
    return false;
  }
  inline bool isOutside(const Vecf<Dim>& pt) const { return isOutside(floatToInt(pt)); }
  // Check if the given cell is free by coordinate
  inline bool isFree(const Veci<Dim>& pn) const {
    if (isOutside(pn))
      return false;
    else
      return isFree(getIndex(pn));
  }
  inline bool isFree(const Vecf<Dim>& pt) const { return isFree(floatToInt(pt)); }
  // Check if the given cell is occupied by coordinate
  inline bool isOccupied(const Veci<Dim>& pn) const {
    if (isOutside(pn))
      return false;
    else
      return isOccupied(getIndex(pn));
  }
  inline bool isOccupied(const Vecf<Dim>& pt) const { return isOccupied(floatToInt(pt)); }
  inline bool isStaticOccupied(const Veci<Dim>& pn) const {
    if (isOutside(pn))
      return false;
    else
      return isStaticOccupied(getIndex(pn));
  }
  inline bool isStaticOccupied(const Vecf<Dim>& pt) const {
    return isStaticOccupied(floatToInt(pt));
  }
  // Check if the given cell is unknown by coordinate
  inline bool isUnknown(const Veci<Dim>& pn) const {
    if (isOutside(pn)) return false;
    return map_[getIndex(pn)] == val_unknown_;
  }
  inline bool isUnknown(const Vecf<Dim>& pt) const { return isUnknown(floatToInt(pt)); }

  // Print basic information about the util
  void info() {
    Vecf<Dim> range = dim_.template cast<decimal_t>() * res_;
    std::cout << "MapUtil Info ========================== " << std::endl;
    std::cout << "   res: [" << res_ << "]" << std::endl;
    std::cout << "   origin: [" << origin_d_.transpose() << "]" << std::endl;
    std::cout << "   range: [" << range.transpose() << "]" << std::endl;
    std::cout << "   dim: [" << dim_.transpose() << "]" << std::endl;
  };

  /** @brief Convert a world-coordinate position to discrete voxel coordinates.
   *  @param pt Position in world coordinates.
   *  @return Integer voxel coordinates.
   */
  inline Veci<Dim> floatToInt(const Vecf<Dim>& pt) const {
    return ((pt - origin_d_) / res_ - Vecf<Dim>::Constant(0.5)).template cast<int>();
  }

  /** @brief Convert discrete voxel coordinates to a world-coordinate position (cell center).
   *  @param pn Integer voxel coordinates.
   *  @return World-coordinate position at the voxel center.
   */
  inline Vecf<Dim> intToFloat(const Veci<Dim>& pn) const {
    return (pn.template cast<decimal_t>() + Vecf<Dim>::Constant(0.5)) * res_ + origin_d_;
  }

  /** @brief Raytrace between two world-coordinate points and return traversed voxels.
   *  @param pt1 Start position.
   *  @param pt2 End position.
   *  @return Ordered vector of unique voxel coordinates along the ray.
   */
  inline vec_Veci<Dim> rayTrace(const Vecf<Dim>& pt1, const Vecf<Dim>& pt2) const {
    Vecf<Dim> diff = pt2 - pt1;
    decimal_t k = 0.1;
    int max_diff = (diff / res_).template lpNorm<Eigen::Infinity>() / k;
    decimal_t s = 1.0 / max_diff;
    Vecf<Dim> step = diff * s;

    vec_Veci<Dim> pns;
    Veci<Dim> prev_pn = Veci<Dim>::Constant(-1);
    for (int n = 1; n < max_diff; n++) {
      Vecf<Dim> pt = pt1 + step * n;
      Veci<Dim> new_pn = floatToInt(pt);
      if (isOutside(new_pn)) break;
      if (new_pn != prev_pn) pns.push_back(new_pn);
      prev_pn = new_pn;
    }
    return pns;
  }

  /** @brief Check if the ray from p1 to p2 is blocked by an occupied voxel.
   *  @param p1 Start position.
   *  @param p2 End position.
   *  @param val Occupancy threshold (voxels with value >= val are considered blocking).
   *  @return True if any voxel along the ray is occupied.
   */
  inline bool isBlocked(const Vecf<Dim>& p1, const Vecf<Dim>& p2, int8_t val = 100) const {
    vec_Veci<Dim> pns = rayTrace(p1, p2);
    for (const auto& pn : pns) {
      if (map_[getIndex(pn)] >= val) return true;
    }
    return false;
  }

  /** @brief Compute the integer bounding box of the vicinity region around a path.
   *  @param path Reference path waypoints.
   *  @param local_box_size Half-extents of the local search box per waypoint.
   *  @param min_point_int Output minimum integer corner of the bounding box.
   *  @param max_point_int Output maximum integer corner of the bounding box.
   */
  void computeVicinityMapInteger(const vec_Vecf<3>& path, const std::vector<float>& local_box_size,
                                 Veci<3>& min_point_int, Veci<3>& max_point_int) const {
    // 1. Compute the global bounding box around the path.
    Vecf<3> min_point_float = Vecf<3>::Constant(std::numeric_limits<float>::max());
    Vecf<3> max_point_float = Vecf<3>::Constant(std::numeric_limits<float>::lowest());

    for (const auto& point : path) {
      // Inflate the local box size (using factor 1.5 as in your example)
      Vecf<3> inflated_local_box_size(1.5 * local_box_size[0], 1.5 * local_box_size[1],
                                      1.5 * local_box_size[2]);
      Vecf<3> local_min = point - inflated_local_box_size;
      Vecf<3> local_max = point + inflated_local_box_size;

      // Update global bounds
      min_point_float = min_point_float.cwiseMin(local_min);
      max_point_float = max_point_float.cwiseMax(local_max);
    }

    // 2. Generate min and max points in integer coordinates
    min_point_int = floatToInt(min_point_float);
    max_point_int = floatToInt(max_point_float);
  }

  // Cloud-getting actually happens here
  template <typename CheckFunc>
  vec_Vecf<Dim> getCloud_(CheckFunc check, const Veci<3>& min_point_int,
                          const Veci<3>& max_point_int) const {
    vec_Vecf<Dim> cloud;

    // Reserve an estimated size (optional, just to reduce reallocations)
    cloud.reserve(static_cast<size_t>((max_point_int - min_point_int).prod()));

    if (Dim == 3) {
#pragma omp parallel
      {
        vec_Vecf<Dim> local_cloud;
// Collapse the three nested loops into one for OpenMP
#pragma omp for collapse(3) nowait
        for (int i = min_point_int(0); i < max_point_int(0); ++i) {
          for (int j = min_point_int(1); j < max_point_int(1); ++j) {
            for (int k = min_point_int(2); k < max_point_int(2); ++k) {
              Veci<3> pti(i, j, k);
              // Use the provided check function
              if (check(pti) && !isOutside(pti)) {
                Vecf<3> ptf = intToFloat(pti);
                local_cloud.push_back(ptf);
              }
            }
          }
        }
#pragma omp critical
        { cloud.insert(cloud.end(), local_cloud.begin(), local_cloud.end()); }
      }
    } else if (Dim == 2) {
#pragma omp parallel
      {
        vec_Vecf<Dim> local_cloud;
#pragma omp for collapse(2) nowait
        for (int i = min_point_int(0); i < max_point_int(0); ++i) {
          for (int j = min_point_int(1); j < max_point_int(1); ++j) {
            Veci<3> pti(i, j, 0);
            if (check(pti) && !isOutside(pti)) {
              Vecf<3> ptf = intToFloat(pti);
              local_cloud.push_back(ptf);
            }
          }
        }
#pragma omp critical
        { cloud.insert(cloud.end(), local_cloud.begin(), local_cloud.end()); }
      }
    }

    return cloud;
  }

  // Cloud-getter
  template <typename CheckFunc>
  vec_Vecf<Dim> getCloud(CheckFunc check, const vec_Vecf<3>& path,
                         const std::vector<float> local_box_size) const {
    vec_Vecf<Dim> cloud;

    // Compute vicinty map integer
    Veci<3> min_point_int, max_point_int;
    computeVicinityMapInteger(path, local_box_size, min_point_int, max_point_int);

    return getCloud_(check, min_point_int, max_point_int);
  }

  // Cloud-getter
  template <typename CheckFunc>
  vec_Vecf<Dim> getCloud(CheckFunc check) const {
    vec_Vecf<Dim> cloud;

    // Get the minimum and maximum points in the current map
    Vecf<3> min_point_float, max_point_float;
    min_point_float(0) = x_min_;
    min_point_float(1) = y_min_;
    min_point_float(2) = z_min_;
    max_point_float(0) = x_max_;
    max_point_float(1) = y_max_;
    max_point_float(2) = z_max_;
    Veci<3> min_point_int = floatToInt(min_point_float);
    Veci<3> max_point_int = floatToInt(max_point_float);

    return getCloud_(check, min_point_int, max_point_int);
  }

  /** @brief Get occupied voxel positions within the vicinity of a path.
   *  @param path Reference path defining the search region.
   *  @param local_box_size Half-extents of the local search box.
   *  @return Vector of occupied voxel positions in world coordinates.
   */
  inline vec_Vecf<Dim> getOccupiedCloud(const vec_Vecf<3>& path,
                                        const std::vector<float> local_box_size) const {
    // Get cloud for occupied cells
    return getCloud([this](const Veci<Dim>& pti) -> bool { return isOccupied(pti); }, path,
                    local_box_size);
  }

  /** @brief Get all occupied voxel positions in the entire map.
   *  @return Vector of occupied voxel positions in world coordinates.
   */
  inline vec_Vecf<Dim> getOccupiedCloud() const {
    // Get cloud for occupied cells
    return getCloud([this](const Veci<Dim>& pti) -> bool { return isOccupied(pti); });
  }

  /** @brief Get free voxel positions within the vicinity of a path.
   *  @param path Reference path defining the search region.
   *  @param local_box_size Half-extents of the local search box.
   *  @return Vector of free voxel positions in world coordinates.
   */
  inline vec_Vecf<Dim> getFreeCloud(const vec_Vecf<3>& path,
                                    const std::vector<float> local_box_size) const {
    // Get cloud for free cells
    return getCloud([this](const Veci<Dim>& pti) -> bool { return isFree(pti); }, path,
                    local_box_size);
  }

  /** @brief Get all free voxel positions in the entire map.
   *  @return Vector of free voxel positions in world coordinates.
   */
  inline vec_Vecf<Dim> getFreeCloud() const {
    // Get cloud for free cells
    return getCloud([this](const Veci<Dim>& pti) -> bool { return isFree(pti); });
  }

  /** @brief Get unknown voxel positions within the vicinity of a path.
   *  @param path Reference path defining the search region.
   *  @param local_box_size Half-extents of the local search box.
   *  @return Vector of unknown voxel positions in world coordinates.
   */
  inline vec_Vecf<Dim> getUnknownCloud(const vec_Vecf<3>& path,
                                       const std::vector<float> local_box_size) const {
    // Get cloud for unknown cells
    return getCloud([this](const Veci<Dim>& pti) -> bool { return isUnknown(pti); }, path,
                    local_box_size);
  }

  /** @brief Get all unknown voxel positions in the entire map.
   *  @return Vector of unknown voxel positions in world coordinates.
   */
  inline vec_Vecf<Dim> getUnknownCloud() const {
    // Get cloud for unknown cells
    return getCloud([this](const Veci<Dim>& pti) -> bool { return isUnknown(pti); });
  }

  /** @brief Get occupied and unknown voxels (treated as occupied) within a path vicinity.
   *  @param path Reference path defining the search region.
   *  @param local_box_size Half-extents of the local search box.
   *  @return Vector of occupied+unknown voxel positions in world coordinates.
   */
  inline vec_Vecf<Dim> getOccupiedCloudWithUnknownAsOccupied(
      const vec_Vecf<3>& path, const std::vector<float> local_box_size) const {
    // Get cloud for occupied cells
    return getCloud(
        [this](const Veci<Dim>& pti) -> bool { return isOccupied(pti) || isUnknown(pti); }, path,
        local_box_size);
  }

  /** @brief Count the number of unknown cells in the map.
   *  @return Number of unknown cells.
   */
  int countUnknownCells() const { return std::count(map_.begin(), map_.end(), val_unknown_); }

  /** @brief Count the number of occupied cells in the map.
   *  @return Number of occupied cells.
   */
  int countOccupiedCells() const { return std::count(map_.begin(), map_.end(), val_occ_); }

  /** @brief Count the number of free cells in the map.
   *  @return Number of free cells.
   */
  int countFreeCells() const { return std::count(map_.begin(), map_.end(), val_free_); }

  /** @brief Get the total number of cells in the voxel map.
   *  @return Total cell count.
   */
  int getTotalNumCells() const { return total_size_; }

  // Map entity
  Tmap map_;

  // ---------------- Dynamic heat-map state ----------------
  std::vector<float> heat_;
  bool dynamic_heat_enabled_{false};
  bool dynamic_as_occupied_current_{
      true};  // Mark current obstacle position as occupied (hard constraint)
  bool dynamic_as_occupied_future_{
      true};            // Mark future reachable positions as occupied (hard constraint)
  float heat_w_{0.0f};  // global planner weight

  // Heat shaping parameters
  float heat_alpha0_{1.0f};
  float heat_alpha1_{2.0f};
  int heat_p_{2};
  int heat_q_{2};
  float heat_tau_ratio_{0.5f};  // tau_w = heat_tau_ratio_ * T_h
  float heat_gamma_{0.0f};      // tube radius growth rate [m/s]
  float heat_Hmax_{10.0f};
  float dyn_base_inflation_m_{0.5f};    // R0 in meters
  float dyn_heat_tube_radius_m_{2.0f};  // Radius of heat corridor [m]
  int heat_num_samples_{15};            // fallback if no prediction samples
  std::vector<vec_Vecf<3>> dyn_pred_samples_;
  std::vector<float> dyn_pred_times_;

  // ---------------- Static obstacle heat parameters ----------------
  bool static_heat_enabled_{false};
  float static_heat_alpha_{2.0f};
  int static_heat_p_{2};
  float static_heat_Hmax_{50.0f};

  // Maximum radius used for precomputing offsets and clamping per-seed function outputs.
  float static_heat_rmax_m_{1.0f};

  // If radius function not set, use this constant radius.
  float static_heat_default_radius_m_{0.5f};

  // Use only boundary occupied voxels as seeds (recommended).
  bool static_heat_boundary_only_{true};

  // Apply halo only on FREE by default; set true to also apply on UNKNOWN.
  bool static_heat_apply_on_unknown_{false};

  // Exclude dynamic-occupied voxels from static heat to avoid double counting.
  bool static_heat_exclude_dynamic_{true};

  // User-provided radius function evaluated at seed voxel center (world coords).
  std::function<float(const Eigen::Vector3f&)> static_heat_radius_fn_;

  // Cached static-heat offsets (recomputed only when resolution/rmax changes).
  struct StaticHeatOff {
    int dx, dy, dz;
    float d_m;
  };
  mutable std::vector<StaticHeatOff> static_heat_off_;
  mutable int static_heat_off_Rcell_{-1};
  mutable float static_heat_off_res_{-1.0f};
  mutable float static_heat_off_rmax_m_{-1.0f};
  mutable std::mutex static_heat_mutex_;

  inline void ensureStaticHeatOffsets(int Rcell) const {
    std::lock_guard<std::mutex> lock(static_heat_mutex_);
    if (static_heat_off_Rcell_ == Rcell && static_heat_off_res_ == (float)res_ &&
        std::fabs(static_heat_off_rmax_m_ - static_heat_rmax_m_) < 1e-6f) {
      return;
    }

    static_heat_off_Rcell_ = Rcell;
    static_heat_off_res_ = (float)res_;
    static_heat_off_rmax_m_ = static_heat_rmax_m_;

    static_heat_off_.clear();
    static_heat_off_.reserve((2 * Rcell + 1) * (2 * Rcell + 1) * (2 * Rcell + 1));

    for (int dx = -Rcell; dx <= Rcell; ++dx) {
      for (int dy = -Rcell; dy <= Rcell; ++dy) {
        for (int dz = -Rcell; dz <= Rcell; ++dz) {
          const float d_m = (float)res_ * std::sqrt(float(dx * dx + dy * dy + dz * dz));
          if (d_m > static_heat_rmax_m_) continue;
          static_heat_off_.push_back({dx, dy, dz, d_m});
        }
      }
    }
  }

  // Soft-cost obstacle mode
  bool use_soft_cost_obstacles_{false};
  float obstacle_soft_cost_{100.0f};

  // Heat cutoff ratio: cells with heat > ratio * Hmax are impassable (0=disabled)
  float heat_cutoff_ratio_{0.5f};

  // ==================== 2D Ground Robot Map ====================

  std::vector<int8_t> map_2d_;         // 2D occupancy grid (dimX * dimY)
  std::vector<float> heat_2d_;         // 2D static heat (dimX * dimY)
  std::vector<float> terrain_cost_;    // 2D terrain gradient cost (dimX * dimY)
  std::vector<float> col_min_height_;  // Per-column min occupied z in world coords
  std::vector<float> col_max_height_;  // Per-column max occupied z in world coords
  bool has_2d_map_{false};

  /** @brief Build a 2D ground robot map from the current 3D occupancy grid.
   *
   *  For each (x,y) column, computes min/max occupied z to determine:
   *  - Hard obstacles: columns with height_span > obstacle_min_height
   *  - Terrain cost: max or avg of height differences with 8-neighbors
   *  - Injects terrain cost into heat_ array at z=0 for A* to use
   *
   *  @param obstacle_min_height Min height span [m] in a column to classify as obstacle.
   *  @param terrain_cost_weight Multiplier for terrain gradient cost injected into heat.
   *  @param cost_mode "max" or "avg" aggregation of neighbor height deltas.
   */
  void buildGroundMap2D(float obstacle_min_height, float terrain_cost_weight,
                        const std::string& cost_mode, bool use_column_any_occupied = true,
                        float column_min_z = 0.0f) {
    if (dim_(0) <= 0 || dim_(1) <= 0 || dim_(2) <= 0) return;

    const int dimX = dim_(0);
    const int dimY = dim_(1);
    const int dimZ = dim_(2);
    const size_t size_2d = static_cast<size_t>(dimX) * dimY;

    // Allocate arrays
    map_2d_.assign(size_2d, val_free_);
    col_min_height_.assign(size_2d, std::numeric_limits<float>::infinity());
    col_max_height_.assign(size_2d, -std::numeric_limits<float>::infinity());
    terrain_cost_.assign(size_2d, 0.0f);

    // Step 1: Column scan — find min/max occupied z per (x,y)
#pragma omp parallel for collapse(2) schedule(static)
    for (int x = 0; x < dimX; ++x) {
      for (int y = 0; y < dimY; ++y) {
        const size_t idx_2d = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y;
        bool has_unknown = false;

        for (int z = 0; z < dimZ; ++z) {
          const size_t idx_3d = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y +
                                static_cast<size_t>(dimX) * static_cast<size_t>(dimY) * z;

          if (map_[idx_3d] == val_occ_) {
            const float world_z = origin_d_(2) + (z + 0.5f) * res_;
            // Skip ground-level voxels below the minimum z threshold
            if (world_z < column_min_z) continue;
            if (world_z < col_min_height_[idx_2d]) col_min_height_[idx_2d] = world_z;
            if (world_z > col_max_height_[idx_2d]) col_max_height_[idx_2d] = world_z;
          } else if (map_[idx_3d] == val_unknown_) {
            has_unknown = true;
          }
        }
      }
    }

    // Step 2: Obstacle classification
#pragma omp parallel for collapse(2) schedule(static)
    for (int x = 0; x < dimX; ++x) {
      for (int y = 0; y < dimY; ++y) {
        const size_t idx_2d = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y;
        if (col_max_height_[idx_2d] > -std::numeric_limits<float>::infinity()) {
          if (use_column_any_occupied) {
            // Any occupied voxel in column → 2D occupied
            map_2d_[idx_2d] = val_occ_;
          } else if (col_min_height_[idx_2d] < std::numeric_limits<float>::infinity()) {
            // Original logic: only mark if height span exceeds threshold
            const float height_span = col_max_height_[idx_2d] - col_min_height_[idx_2d];
            if (height_span > obstacle_min_height) {
              map_2d_[idx_2d] = val_occ_;
            }
          }
        }
      }
    }

    // Step 3: Terrain cost computation for free cells
    const bool use_max = (cost_mode == "max");
#pragma omp parallel for collapse(2) schedule(static)
    for (int x = 0; x < dimX; ++x) {
      for (int y = 0; y < dimY; ++y) {
        const size_t idx_2d = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y;
        if (map_2d_[idx_2d] != val_free_) continue;

        // Skip cells with no height data
        const float h_center = col_max_height_[idx_2d];
        if (h_center == -std::numeric_limits<float>::infinity()) continue;

        float sum_delta = 0.0f;
        float max_delta = 0.0f;
        int count = 0;

        for (int dx = -1; dx <= 1; ++dx) {
          for (int dy = -1; dy <= 1; ++dy) {
            if (dx == 0 && dy == 0) continue;
            const int nx = x + dx;
            const int ny = y + dy;
            if (nx < 0 || nx >= dimX || ny < 0 || ny >= dimY) continue;

            const size_t n_idx = static_cast<size_t>(nx) + static_cast<size_t>(dimX) * ny;
            const float h_neigh = col_max_height_[n_idx];
            if (h_neigh == -std::numeric_limits<float>::infinity()) continue;

            const float delta = std::abs(h_center - h_neigh);
            sum_delta += delta;
            if (delta > max_delta) max_delta = delta;
            ++count;
          }
        }

        if (count > 0) {
          terrain_cost_[idx_2d] = use_max ? max_delta : (sum_delta / static_cast<float>(count));
        }
      }
    }

    // Step 4: Inject terrain cost into heat_ array at z=0 slice
    if (terrain_cost_weight > 0.0f && !heat_.empty() && dimZ > 0) {
      for (int x = 0; x < dimX; ++x) {
        for (int y = 0; y < dimY; ++y) {
          const size_t idx_2d = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y;
          if (terrain_cost_[idx_2d] > 0.0f) {
            // z=0 slice in the 3D heat array
            const size_t idx_3d = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y;
            if (idx_3d < heat_.size()) {
              heat_[idx_3d] += terrain_cost_weight * terrain_cost_[idx_2d];
            }
          }
        }
      }
    }

    // Step 5: Compute 2D static heat halos around occupied cells
    heat_2d_.assign(size_2d, 0.0f);
    if (static_heat_enabled_) {
      const int Rcell = static_cast<int>(std::ceil(static_heat_rmax_m_ / res_));
      // Precompute 2D offsets within radius
      struct Off2D {
        int dx, dy;
        float d_m;
      };
      std::vector<Off2D> offsets_2d;
      for (int dx = -Rcell; dx <= Rcell; ++dx) {
        for (int dy = -Rcell; dy <= Rcell; ++dy) {
          float d_m = res_ * std::sqrt(float(dx * dx + dy * dy));
          if (d_m <= static_heat_rmax_m_) offsets_2d.push_back({dx, dy, d_m});
        }
      }

      for (int x0 = 0; x0 < dimX; ++x0) {
        for (int y0 = 0; y0 < dimY; ++y0) {
          const size_t src_idx = static_cast<size_t>(x0) + static_cast<size_t>(dimX) * y0;
          if (map_2d_[src_idx] != val_occ_) continue;

          for (const auto& o : offsets_2d) {
            const int x = x0 + o.dx, y = y0 + o.dy;
            if (x < 0 || x >= dimX || y < 0 || y >= dimY) continue;
            const size_t idx = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y;

            // Distance-decay heat
            const float u = std::min(1.0f, std::max(0.0f, o.d_m / static_heat_rmax_m_));
            const float base = 1.0f - u;
            float pw = (static_heat_p_ == 2)   ? base * base
                       : (static_heat_p_ == 3) ? base * base * base
                                               : std::pow(base, float(static_heat_p_));
            float w = std::min(static_heat_alpha_ * pw, static_heat_Hmax_);
            if (w > heat_2d_[idx]) heat_2d_[idx] = w;
          }
        }
      }
    }

    // Step 6: Project 3D dynamic heat into 2D (max over z columns)
    if (dynamic_heat_enabled_ && !heat_.empty()) {
      for (int x = 0; x < dimX; ++x) {
        for (int y = 0; y < dimY; ++y) {
          float max_h = 0.0f;
          for (int z = 0; z < dimZ; ++z) {
            const size_t idx_3d =
                static_cast<size_t>(x) + dimX * (static_cast<size_t>(y) + dimY * z);
            if (idx_3d < heat_.size()) max_h = std::max(max_h, heat_[idx_3d]);
          }
          const size_t idx_2d = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y;
          if (max_h > heat_2d_[idx_2d]) heat_2d_[idx_2d] = max_h;
        }
      }
    }

    has_2d_map_ = true;
  }

  /** @brief Check if a 2D ground map has been built. */
  bool has2DMap() const { return has_2d_map_; }

  /** @brief Get 2D occupancy at grid coordinates. */
  int8_t get2DOccupancy(int x, int y) const {
    if (!has_2d_map_) return val_free_;
    const int dimX = dim_(0);
    const int dimY = dim_(1);
    if (x < 0 || x >= dimX || y < 0 || y >= dimY) return val_occ_;
    return map_2d_[static_cast<size_t>(x) + static_cast<size_t>(dimX) * y];
  }

  /** @brief Get terrain gradient cost at grid coordinates. */
  float getTerrainCost(int x, int y) const {
    if (!has_2d_map_) return 0.0f;
    const int dimX = dim_(0);
    const int dimY = dim_(1);
    if (x < 0 || x >= dimX || y < 0 || y >= dimY) return 0.0f;
    return terrain_cost_[static_cast<size_t>(x) + static_cast<size_t>(dimX) * y];
  }

  /** @brief Get terrain height (max occupied z) at grid coordinates.
   *  Returns 0 if no occupied voxels found in that column. */
  float getTerrainHeight(int x, int y) const {
    if (!has_2d_map_) return 0.0f;
    const int dimX = dim_(0);
    const int dimY = dim_(1);
    if (x < 0 || x >= dimX || y < 0 || y >= dimY) return 0.0f;
    const float h = col_max_height_[static_cast<size_t>(x) + static_cast<size_t>(dimX) * y];
    return (h == -std::numeric_limits<float>::infinity()) ? 0.0f : h;
  }

  /** @brief Get terrain height at world coordinates with nearest-cell lookup.
   *  @param wx World x coordinate.
   *  @param wy World y coordinate.
   *  @return Terrain height at that location, or 0 if unavailable.
   */
  float getTerrainHeightWorld(float wx, float wy) const {
    if (!has_2d_map_) return 0.0f;
    const int ix = static_cast<int>(std::floor((wx - origin_d_(0)) / res_ - 0.5f));
    const int iy = static_cast<int>(std::floor((wy - origin_d_(1)) / res_ - 0.5f));
    return getTerrainHeight(ix, iy);
  }

  /** @brief Get raw pointer to the 2D occupancy map data for GraphSearch construction. */
  const int8_t* get2DMapData() const { return has_2d_map_ ? map_2d_.data() : nullptr; }

  /** @brief Get raw pointer to the 2D heat data for GraphSearch. */
  const float* get2DHeatData() const { return has_2d_map_ ? heat_2d_.data() : nullptr; }

  /** @brief Free a cell and its surroundings in the 2D map. */
  void free2DCell(int cx, int cy, float radius_m) {
    if (!has_2d_map_) return;
    const int dimX = dim_(0);
    const int dimY = dim_(1);
    const int r = static_cast<int>(std::ceil(radius_m / res_));
    for (int dx = -r; dx <= r; ++dx) {
      for (int dy = -r; dy <= r; ++dy) {
        const int nx = cx + dx;
        const int ny = cy + dy;
        if (nx >= 0 && nx < dimX && ny >= 0 && ny < dimY) {
          map_2d_[static_cast<size_t>(nx) + static_cast<size_t>(dimX) * ny] = val_free_;
        }
      }
    }
  }

  /** @brief Get 2D map dimensions. */
  void get2DDimensions(int& dimX, int& dimY) const {
    dimX = dim_(0);
    dimY = dim_(1);
  }

  /** @brief Build 2D map entirely from ESDF grid. Trusts ESDF for both occupancy and heat.
   *
   *  - Occupied: ESDF distance = 0 (on obstacle in mapper's grid)
   *  - Heat: linear ramp from h_max at d=0 → 0 at d=d_safe
   *  - No point cloud processing or inflation involved.
   *
   *  @param esdf   The 2D ESDF grid from the mapper.
   *  @param d_safe Distance at which heat reaches 0.
   *  @param h_max  Maximum heat value (at obstacle surface).
   */
  /** @brief Project 3D dynamic heat (heat_) into 2D heat_2d_ as max-over-z.
   *
   *  Both ESDF/Occ2D 2D builders compute heat_2d_ from static distance fields
   *  only; they would otherwise discard the dynamic-obstacle tube populated in
   *  heat_ during updateMap. Calling this after their static fill merges the
   *  dynamic tube in via per-column max, so 2D A* (and ground_2d_heat viz)
   *  account for predicted-trajectory heat.
   */
  void mergeDynamicHeatInto2D() {
    if (!dynamic_heat_enabled_ || heat_.empty() || heat_2d_.empty()) return;
    const int dimX = dim_(0);
    const int dimY = dim_(1);
    const int dimZ = dim_(2);
    if (dimX <= 0 || dimY <= 0 || dimZ <= 0) return;
    for (int x = 0; x < dimX; ++x) {
      for (int y = 0; y < dimY; ++y) {
        float max_h = 0.0f;
        for (int z = 0; z < dimZ; ++z) {
          const size_t idx_3d = static_cast<size_t>(x) +
                                static_cast<size_t>(dimX) *
                                    (static_cast<size_t>(y) + static_cast<size_t>(dimY) * z);
          if (idx_3d < heat_.size()) max_h = std::max(max_h, heat_[idx_3d]);
        }
        const size_t idx_2d = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y;
        if (max_h > heat_2d_[idx_2d]) heat_2d_[idx_2d] = max_h;
      }
    }
  }

  void buildMap2DFromEsdf(const EsdfGrid2D& esdf, double d_safe, double h_max) {
    const int dimX = dim_(0);
    const int dimY = dim_(1);
    const size_t n2d = static_cast<size_t>(dimX) * dimY;

    map_2d_.assign(n2d, val_free_);
    heat_2d_.assign(n2d, 0.0f);

    // A MIGHTY cell of width res_ centered at (wx, wy) intersects an obstacle
    // if the ESDF distance anywhere inside the cell is <= 0. Querying only the
    // center misses obstacles that clip the cell's corners. The half-diagonal
    // of the cell (res_*sqrt(2)/2) is the worst-case distance from the center
    // to any point inside, so anything with center distance <= half_diag could
    // be occupied somewhere within the cell.
    const double half_diag = res_ * 0.5 * std::sqrt(2.0);

    for (int y = 0; y < dimY; ++y) {
      for (int x = 0; x < dimX; ++x) {
        const double wx = origin_d_(0) + (x + 0.5) * res_;
        const double wy = origin_d_(1) + (y + 0.5) * res_;
        const size_t idx = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y;

        if (esdf.isInBounds(wx, wy)) {
          const double d = esdf.queryDistance(wx, wy);
          if (d <= half_diag) {
            map_2d_[idx] = val_occ_;
            heat_2d_[idx] = static_cast<float>(h_max);
          } else if (d < d_safe) {
            heat_2d_[idx] = static_cast<float>(h_max * (1.0 - d / d_safe));
          }
        }
        // Outside ESDF bounds: leave as free (initialized above). Marking
        // OOB as occupied creates phantom walls when MIGHTY's robot-centered
        // window extends past the mapper's fixed-origin coverage.
      }
    }
    mergeDynamicHeatInto2D();
    has_2d_map_ = true;
  }

  /** @brief Build 2D map from binary occupancy grid with BFS distance-based heat.
   *
   *  Occupancy comes directly from the mapper's occ_2d_topic (binary 0/100).
   *  Heat is computed via BFS distance transform: h_max at obstacle → 0 at d_safe.
   *
   *  @param occ    Binary 2D occupancy grid from mapper.
   *  @param d_safe Distance at which heat reaches 0 [m].
   *  @param h_max  Maximum heat value at occupied cells.
   */
  void buildMap2DFromOcc2D(const OccGrid2D& occ, double d_safe, double h_max) {
    const int dimX = dim_(0);
    const int dimY = dim_(1);
    const size_t n2d = static_cast<size_t>(dimX) * dimY;

    map_2d_.assign(n2d, val_free_);
    heat_2d_.assign(n2d, 0.0f);

    // Compute distance field from the binary occupancy grid (truncated at d_safe)
    std::vector<float> dist = occ.computeDistanceField(d_safe);

    const double inv_occ_res = occ.invResolution();
    const int occ_w = occ.width();
    const int occ_h = occ.height();
    const auto& occ_data = occ.occupiedData();

    for (int y = 0; y < dimY; ++y) {
      for (int x = 0; x < dimX; ++x) {
        // Convert MIGHTY grid cell to world coordinates (cell extent is [wx-h, wx+h])
        const double wx = origin_d_(0) + (x + 0.5) * res_;
        const double wy = origin_d_(1) + (y + 0.5) * res_;
        const double half = 0.5 * res_;
        const size_t idx = static_cast<size_t>(x) + static_cast<size_t>(dimX) * y;

        // Range of underlying occ-grid cells overlapping this MIGHTY cell.
        int ox_min = static_cast<int>(std::floor((wx - half - occ.originX()) * inv_occ_res));
        int ox_max = static_cast<int>(std::floor((wx + half - occ.originX()) * inv_occ_res));
        int oy_min = static_cast<int>(std::floor((wy - half - occ.originY()) * inv_occ_res));
        int oy_max = static_cast<int>(std::floor((wy + half - occ.originY()) * inv_occ_res));

        // Clip to source grid. Cells with no overlap at all stay free
        // (don't ring the map with phantom walls when MIGHTY's window
        // extends past the mapper's fixed-origin coverage).
        ox_min = std::max(ox_min, 0);
        oy_min = std::max(oy_min, 0);
        ox_max = std::min(ox_max, occ_w - 1);
        oy_max = std::min(oy_max, occ_h - 1);
        if (ox_min > ox_max || oy_min > oy_max) {
          continue;
        }

        // Mark occupied if ANY underlying source cell is occupied. Otherwise
        // pull heat from the nearest underlying source cell (smallest d → max heat).
        // NOTE: occ_data only flags true-occupied cells (>=100). Unknown cells
        // (-1) are implicitly traversable for HGP A* — this gives us the
        // "plan through unknown" semantics that frontier exploration needs.
        bool any_occ = false;
        float min_d = static_cast<float>(d_safe);
        for (int oy = oy_min; oy <= oy_max && !any_occ; ++oy) {
          for (int ox = ox_min; ox <= ox_max; ++ox) {
            const size_t oidx = static_cast<size_t>(oy) * occ_w + ox;
            if (occ_data[oidx]) {
              any_occ = true;
              break;
            }
            if (dist[oidx] < min_d) min_d = dist[oidx];
          }
        }

        if (any_occ) {
          map_2d_[idx] = val_occ_;
          heat_2d_[idx] = static_cast<float>(h_max);
        } else if (min_d < static_cast<float>(d_safe)) {
          heat_2d_[idx] = static_cast<float>(h_max * (1.0 - min_d / d_safe));
        }
      }
    }
    mergeDynamicHeatInto2D();
    has_2d_map_ = true;
  }

 protected:
  // Resolution
  decimal_t res_;
  // Total size of the map
  int total_size_ = 0;
  // Inflation
  float inflation_;
  // Origin, float type
  Vecf<Dim> origin_d_;
  // Center, float type
  Vecf<Dim> center_map_;
  // Dimension, int type
  Veci<Dim> dim_, prev_dim_;
  // Cached stride for 3D indexing (dim[0] * dim[1])
  int dim_xy_ = 0;
  // Map values
  float x_map_min_, x_map_max_, y_map_min_, y_map_max_, z_map_min_, z_map_max_;
  float x_min_, x_max_, y_min_, y_max_, z_min_, z_max_;
  // Obstacle maximum velocity
  float obst_max_vel_;
  // Cells size
  int cells_x_, cells_y_, cells_z_;
  // Assume occupied cell has value 100
  int8_t val_occ_ = 100;
  // Assume free cell has value 0
  int8_t val_free_ = 0;
  // Assume unknown cell has value -1
  int8_t val_unknown_ = -1;

  // Flags
  bool map_initialized_ = false;

  // small map buffer
  Vec3f min_point_;
  Vec3f max_point_;
};

typedef MapUtil<3> VoxelMapUtil;
}  // namespace mighty