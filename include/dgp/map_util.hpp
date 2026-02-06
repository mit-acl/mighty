/**
 * @file map_util.h
 * @brief MapUtil classes
 */
#ifndef DGP_MAP_UTIL_H
#define DGP_MAP_UTIL_H

#include <iostream>
#include "dgp/data_type.hpp"
#include <mighty/mighty_type.hpp>
#include <pcl/kdtree/kdtree_flann.h>
#include "timer.hpp"
#include <omp.h>

namespace mighty
{

  // The type of map data Tmap is defined as a 1D array
  // Using int8_t saves 75% memory compared to int (only need 3 values: -1, 0, 100)
  using Tmap = std::vector<int8_t>;
  typedef timer::Timer MyTimer;

  /**
   * @brief The map util class for collision checking
   * @param Dim is the dimension of the workspace
   */
  template <int Dim>
  class MapUtil
  {
  public:
    // Constructor
    MapUtil(float res, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, float inflation, float obst_max_vel = 1.0f)
    {

      /* --------- Initialize parameters --------- */
      setInflation(inflation);                                                                                                                               // Set inflation
      setResolution(res);                                                                                                                                    // Set the resolution
      setMapSize(x_min, x_max, y_min, y_max, z_min, z_max);                                                                                                  // Set the cells and z_boundaries
      obst_max_vel_ = obst_max_vel;
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
        obst_max_vel_(other.obst_max_vel_),
        res_(other.res_),
        total_size_(other.total_size_),
        inflation_(other.inflation_),
        origin_d_(other.origin_d_),
        center_map_(other.center_map_),
        dim_(other.dim_),
        prev_dim_(other.prev_dim_),
        dim_xy_(other.dim_xy_),
        x_map_min_(other.x_map_min_), x_map_max_(other.x_map_max_),
        y_map_min_(other.y_map_min_), y_map_max_(other.y_map_max_),
        z_map_min_(other.z_map_min_), z_map_max_(other.z_map_max_),
        x_min_(other.x_min_), x_max_(other.x_max_),
        y_min_(other.y_min_), y_max_(other.y_max_),
        z_min_(other.z_min_), z_max_(other.z_max_),
        cells_x_(other.cells_x_), cells_y_(other.cells_y_), cells_z_(other.cells_z_),
        val_occ_(other.val_occ_), val_free_(other.val_free_), val_unknown_(other.val_unknown_),
        map_initialized_(other.map_initialized_),
        min_point_(other.min_point_),
        max_point_(other.max_point_)
    {
      // Mutex is default-constructed
    }

    // Heat map configuration
    void setObstMaxVelocity(float v) { obst_max_vel_ = v; }
    void setDynamicHeatEnabled(bool e) { dynamic_heat_enabled_ = e; }
    void setStaticHeatEnabled(bool e) { static_heat_enabled_ = e; }
    void setHeatWeight(float w) { heat_w_ = w; }
    bool dynamicHeatEnabled() const { return dynamic_heat_enabled_; }
    bool staticHeatEnabled() const { return static_heat_enabled_; }
    float getHeatWeight() const { return heat_w_; }

    void setDynamicHeatParams(float alpha0, float alpha1, int p, int q,
                              float tau_ratio, float gamma, float Hmax, float base_infl)
    {
      heat_alpha0_ = alpha0;
      heat_alpha1_ = alpha1;
      heat_p_ = p;
      heat_q_ = q;
      heat_tau_ratio_ = tau_ratio;
      heat_gamma_ = gamma;
      heat_Hmax_ = Hmax;
      dyn_base_inflation_m_ = base_infl;
    }

    void setDynamicPredictedSamples(const std::vector<vec_Vecf<3>> &samples,
                                    const std::vector<float> &times)
    {
      dyn_pred_samples_ = samples;
      dyn_pred_times_ = times;
    }

    void setStaticHeatParams(float alpha, int p, float Hmax, float rmax,
                            bool boundary_only = true, bool on_unknown = false,
                            bool exclude_dyn = true)
    {
      static_heat_alpha_ = alpha;
      static_heat_p_ = p;
      static_heat_Hmax_ = Hmax;
      static_heat_rmax_m_ = rmax;
      static_heat_boundary_only_ = boundary_only;
      static_heat_apply_on_unknown_ = on_unknown;
      static_heat_exclude_dynamic_ = exclude_dyn;
    }

    void setStaticHeatRadiusFunction(
        const std::function<float(const Eigen::Vector3f &)> &fn,
        float default_radius_m)
    {
      static_heat_radius_fn_ = fn;
      static_heat_default_radius_m_ = default_radius_m;
    }

    // Query heat at voxel coordinates
    float getHeat(int x, int y, int z) const
    {
      if (heat_.empty() || x < 0 || x >= dim_(0) ||
          y < 0 || y >= dim_(1) || z < 0 || z >= dim_(2))
        return 0.0f;
      const int idx = x + y * dim_(0) + z * dim_xy_;
      return (idx >= 0 && idx < (int)heat_.size()) ? heat_[idx] : 0.0f;
    }

    // Get heat cloud for visualization
    vec_Vecf<3> getHeatCloud(float threshold = 0.01f) const;
    std::vector<float> getHeatValues() const { return heat_; }
    float getMaxHeat() const;

    // Destructor
    ~MapUtil()
    {
      // Clear the map
      map_.clear();
    }

    // assume Vec3f is Eigen::Vector3f and Vec3i is Eigen::Vector3i
    void readMap(
        const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& cloud,
        int cells_x, int cells_y, int cells_z,
        const Vec3f &center_map,
        double z_ground,
        double z_max,
        double inflation,
        const vec_Vecf<3> &obst_pos = vec_Vecf<3>(),
        const vec_Vecf<3> &obst_bbox = vec_Vecf<3>(),
        double traj_max_time = 0.0)
    {
      // 1) Compute X/Y dims with inflation pad
      int pad = int(std::ceil(5.0 * inflation / res_));
      int dimX = cells_x + pad, dimY = cells_y + pad, dimZ = cells_z;

      // 2) Compute how many cells below/above center we keep,
      //    strictly within [z_ground, z_max]
      int halfZ = dimZ / 2;
      int down = halfZ, up = halfZ;
      // world coords of bottom slice:
      float bot = center_map.z() - halfZ * res_;
      if (bot < z_ground)
        down = std::max(int(std::floor((center_map.z() - z_ground) / res_)), 0);
      // top slice:
      float top = center_map.z() + halfZ * res_;
      if (top > z_max)
        up = std::max(int(std::floor((z_max - center_map.z()) / res_)), 1);
      dimZ = down + up;

      // 3) Compute origin (global coords of cell (0,0,0)) and clamp it
      Vec3f origin;
      origin.x() = center_map.x() - (dimX * res_) / 2.0f;
      origin.y() = center_map.y() - (dimY * res_) / 2.0f;
      origin.z() = center_map.z() - down * res_;
      // ensure origin.z >= z_ground and origin.z+dimZ*res <= z_max
      origin.z() = std::clamp(origin.z(),
                              z_ground,
                              z_max - dimZ * res_);

      // 4) Allocate map and fill as unknown
      size_t total = size_t(dimX) * dimY * dimZ;
      map_.assign(total, val_unknown_);

      // 5) Precompute inflation offsets
      int m = int(std::floor(inflation / res_));
      std::vector<Vec3i> offsets;
      offsets.reserve((2 * m + 1) * (2 * m + 1) * (2 * m + 1));
      for (int dx = -m; dx <= m; ++dx)
        for (int dy = -m; dy <= m; ++dy)
          for (int dz = -m; dz <= m; ++dz)
            offsets.emplace_back(dx, dy, dz);

      // 6) Helpers for clamping & indexing
      auto clamp_idx = [&](int v, int M)
      { return std::clamp(v, 0, M - 1); };
      auto idx3 = [&](int x, int y, int z)
      {
        return size_t(x) + size_t(dimX) * y + size_t(dimX) * size_t(dimY) * z;
      };

      // 7) Rasterize & inflate, skipping points outside z-bounds
      #pragma omp parallel for schedule(dynamic)
      for (size_t i = 0; i < cloud->points.size(); ++i)
      {
        const auto &P = cloud->points[i];
        if (P.z < z_ground || P.z > z_max)
          continue;
        int xi = clamp_idx(int(std::floor((P.x - origin.x()) / res_)), dimX);
        int yi = clamp_idx(int(std::floor((P.y - origin.y()) / res_)), dimY);
        int zi = clamp_idx(int(std::floor((P.z - origin.z()) / res_)), dimZ);

        // mark occupied
        map_[idx3(xi, yi, zi)] = val_occ_;
        // inflate neighborhood
        for (auto &off : offsets)
        {
          int x2 = xi + off.x(), y2 = yi + off.y(), z2 = zi + off.z();
          if (x2 < 0 || x2 >= dimX || y2 < 0 || y2 >= dimY || z2 < 0 || z2 >= dimZ)
            continue;
          map_[idx3(x2, y2, z2)] = val_occ_;
        }
      }

      // 8) Update metadata
      dim_ = Veci<3>(dimX, dimY, dimZ);
      dim_xy_ = dimX * dimY; // Cache the stride for fast indexing
      total_size_ = total;
      origin_d_ = origin;
      center_map_ = center_map;

      // 9) Compute heat maps if enabled
      if (dynamic_heat_enabled_ || static_heat_enabled_)
      {
        // Initialize heat vector to same size as map
        heat_.assign(total, 0.0f);

        // Compute dynamic heat from moving obstacles
        if (dynamic_heat_enabled_)
        {
          computeDynamicHeat(obst_pos, obst_bbox, traj_max_time);
        }

        // Compute static heat (boundary halo)
        if (static_heat_enabled_)
        {
          computeStaticHeat();
        }
      }
      else
      {
        // Clear heat if disabled
        heat_.clear();
      }
    }


    // Pre-compute inflation
    // Precompute offsets for inflation
    vec_Veci<3> computeInflationOffsets(const Veci<3> &inflation_cells)
    {
      vec_Veci<3> offsets;

      // include diagonal offsets
      for (int dx = -inflation_cells(0); dx <= inflation_cells(0); ++dx)
      {
        for (int dy = -inflation_cells(1); dy <= inflation_cells(1); ++dy)
        {
          for (int dz = -inflation_cells(2); dz <= inflation_cells(2); ++dz)
          {
            offsets.push_back(Veci<3>(dx, dy, dz));
          }
        }
      }

      return offsets;
    }

    Veci<3> indexToVeci3(int index)
    {
      Veci<3> position;
      position[0] = index % dim_(0);
      position[1] = (index / dim_(0)) % dim_(1);
      position[2] = index / (dim_(0) * dim_(1));
      return position;
    }

    void setCellSize(int cells_x, int cells_y, int cells_z)
    {
      // Set cells
      cells_x_ = cells_x;
      cells_y_ = cells_y;
      cells_z_ = cells_z;
    }

    void setMapSize(float x_min, float x_max, float y_min, float y_max, float z_min, float z_max)
    {
      // Set map boundaries
      x_map_min_ = x_min;
      x_map_max_ = x_max;
      y_map_min_ = y_min;
      y_map_max_ = y_max;
      z_map_min_ = z_min;
      z_map_max_ = z_max;
    }

    /**
     * @brief  Find a free point in the map that is closest to the given point
     * @param  vec_Vecf<3> point : The given point
     * @param  vec_Vecf<3> free_point : The free point that is closest to the given point
     * @return void
     */
    void findClosestFreePoint(const Vec3f &point, Vec3f &closest_free_point)
    {

      // Initialize the closest free point
      closest_free_point = point;

      // Check if the map is initialized
      if (!map_initialized_)
      {
        std::cout << "Map is not initialized" << std::endl;
        return;
      }

      // Get the position of the point in int
      Veci<3> point_int = floatToInt(point);

      // Get the index
      int index = getIndex(point_int);

      if (index >= 0 && index < total_size_)
      {
        // Check if the point is free
        if (map_[index] == val_free_)
        {
          closest_free_point = point;
          return;
        }

        // Get the neighboring indices - preallocate to avoid repeated allocations
        std::vector<int> neighbor_indices;
        neighbor_indices.reserve(1000); // Pre-allocate for typical case

        // Increase the radius until a free point is found
        for (float radius = 1.0; radius < 5.0; radius += 0.5) // TODO: expose the radius as a parameter
        {
          neighbor_indices.clear();
          getNeighborIndicesSphere(point_int, neighbor_indices, radius);

          // Find the closest free point
          float min_dist = std::numeric_limits<float>::max();
          for (int neighbor_index : neighbor_indices)
          {
            if (neighbor_index >= 0 && neighbor_index < total_size_)
            {
              if (map_[neighbor_index] == val_free_)
              {
                Veci<3> neighbor_int = indexToVeci3(neighbor_index);
                Vec3f neighbor = intToFloat(neighbor_int);
                float dist = (neighbor - point).norm();
                if (dist < min_dist)
                {
                  min_dist = dist;
                  closest_free_point = neighbor;
                }
              }
            }
          }

          // Check if a free point is found
          if (min_dist < std::numeric_limits<float>::max())
          {
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
    void getNeighborIndices(const Veci<3> &point_int, std::vector<int> &neighbor_indices, float radius)
    {
      // Get the radius in int
      float radius_int = radius / res_;
      Veci<3> radius_int_vec(radius_int, radius_int, radius_int);

      // Get the min and max positions
      Veci<3> min_pos = point_int - radius_int_vec;
      Veci<3> max_pos = point_int + radius_int_vec;

      // Iterate over the neighbors
      for (int x = min_pos[0]; x <= max_pos[0]; ++x)
      {
        for (int y = min_pos[1]; y <= max_pos[1]; ++y)
        {
          for (int z = min_pos[2]; z <= max_pos[2]; ++z)
          {

            // Check if the neighbor is inside the map
            if (x >= 0 && x < dim_(0) && y >= 0 && y < dim_(1) && z >= 0 && z < dim_(2))
            {

              Veci<3> neighbor_int(x, y, z);
              int index = getIndex(neighbor_int);
              if (index >= 0 && index < total_size_)
              {
                neighbor_indices.push_back(index);
              }
            }
          }
        }
      }
    }

    /**
     * @brief Get indices of neighbors within spherical radius (optimized)
     * @param Veci<3> point_int : Center point
     * @param std::vector<int>& neighbor_indices : Output indices
     * @param float radius : Spherical radius
     * @return void
     */
    void getNeighborIndicesSphere(const Veci<3> &point_int, std::vector<int> &neighbor_indices, float radius)
    {
      // Get the radius in int
      float radius_int = radius / res_;
      int r_int = static_cast<int>(std::ceil(radius_int));
      float r_sq = radius_int * radius_int;

      // Get the bounding box
      int x_min = std::max(0, point_int[0] - r_int);
      int x_max = std::min(dim_(0) - 1, point_int[0] + r_int);
      int y_min = std::max(0, point_int[1] - r_int);
      int y_max = std::min(dim_(1) - 1, point_int[1] + r_int);
      int z_min = std::max(0, point_int[2] - r_int);
      int z_max = std::min(dim_(2) - 1, point_int[2] + r_int);

      // Iterate over the bounding box and check spherical distance
      for (int x = x_min; x <= x_max; ++x)
      {
        int dx = x - point_int[0];
        int dx_sq = dx * dx;
        for (int y = y_min; y <= y_max; ++y)
        {
          int dy = y - point_int[1];
          int dxy_sq = dx_sq + dy * dy;
          for (int z = z_min; z <= z_max; ++z)
          {
            int dz = z - point_int[2];
            float dist_sq = dxy_sq + dz * dz;

            // Only include points within spherical radius
            if (dist_sq <= r_sq)
            {
              Veci<3> neighbor_int(x, y, z);
              int index = getIndex(neighbor_int);
              if (index >= 0 && index < total_size_)
              {
                neighbor_indices.push_back(index);
              }
            }
          }
        }
      }
    }

    // Check if the given point has any occupied neighbors.
    // Returns true if at least one neighboring cell is non-free.
    inline bool checkIfPointHasNonFreeNeighbour(const Veci<Dim> &pt) const
    {
      if constexpr (Dim == 2)
      {
        for (int dx = -1; dx <= 1; ++dx)
        {
          for (int dy = -1; dy <= 1; ++dy)
          {
            // Skip the center point
            if (dx == 0 && dy == 0)
              continue;
            Veci<2> neighbor = pt;
            neighbor(0) += dx;
            neighbor(1) += dy;
            // Check if the neighbor is within the map and non-free.
            if (!isOutside(neighbor) && !isFree(neighbor))
              return true;
          }
        }
      }
      else if constexpr (Dim == 3)
      {
        for (int dx = -1; dx <= 1; ++dx)
        {
          for (int dy = -1; dy <= 1; ++dy)
          {
            for (int dz = -1; dz <= 1; ++dz)
            {
              // Skip the center point
              if (dx == 0 && dy == 0 && dz == 0)
                continue;
              Veci<3> neighbor = pt;
              neighbor(0) += dx;
              neighbor(1) += dy;
              neighbor(2) += dz;
              // Check if the neighbor is within the map and non-free.
              if (!isOutside(neighbor) && !isFree(neighbor))
                return true;
            }
          }
        }
      }
      return false;
    }

    void setInflation(float inflation)
    {
      inflation_ = inflation;
    }

    void setResolution(float res)
    {
      res_ = res;
    }

    // Check if a point is inside the box formed by 8 points (convex hull of 8 points)
    bool isPointInBox(const Vec3f &p, const std::vector<Vec3f> &vertices)
    {
      // Assuming vertices contains 8 points representing the corners of the box in 3D space.
      // The box is axis-aligned and the vertices are provided in any order.

      // Find the minimum and maximum x, y, and z coordinates among the vertices.
      Vec3f min = vertices[0];
      Vec3f max = vertices[0];

      for (const auto &vertex : vertices)
      {
        min.x() = std::min(min.x(), vertex.x());
        min.y() = std::min(min.y(), vertex.y());
        min.z() = std::min(min.z(), vertex.z());

        max.x() = std::max(max.x(), vertex.x());
        max.y() = std::max(max.y(), vertex.y());
        max.z() = std::max(max.z(), vertex.z());
      }

      // Add buffer to the min and max coordinates to give some margin around the box.
      // TODO: expose this buffer as a parameter.
      min -= Vec3f(0.5, 0.5, 0.5); // min -= Vec3f
      max += Vec3f(0.5, 0.5, 0.5); // max += Vec3f

      // Check if the point lies within the bounds defined by the min and max coordinates.
      return (p.x() >= min.x() && p.x() <= max.x()) &&
             (p.y() >= min.y() && p.y() <= max.y()) &&
             (p.z() >= min.z() && p.z() <= max.z());
    }

    // Get resolution
    decimal_t getRes()
    {
      return res_;
    }
    // Get dimensions
    Veci<Dim> getDim()
    {
      return dim_;
    }
    // Get origin
    Vecf<Dim> getOrigin()
    {
      return origin_d_;
    }
    // Get index of a cell
    inline int getIndex(const Veci<Dim> &pn) const
    {
      return Dim == 2 ? pn(0) + dim_(0) * pn(1) : pn(0) + dim_(0) * pn(1) + dim_xy_ * pn(2);
    }
    // Get index of a cell in old map
    inline int getOldIndex(const Veci<Dim> &pn) const
    {
      return Dim == 2 ? pn(0) + prev_dim_(0) * pn(1) : pn(0) + prev_dim_(0) * pn(1) + prev_dim_(0) * prev_dim_(1) * pn(2);
    }

    ///
    Veci<Dim> getVoxelPos(int idx)
    {
      Veci<Dim> pn;
      if (Dim == 2)
      {
        pn(0) = idx % dim_(0);
        pn(1) = idx / dim_(0);
      }
      else
      {
        pn(0) = idx % dim_(0);
        pn(1) = (idx / dim_(0)) % dim_(1);
        pn(2) = idx / (dim_(0) * dim_(1));
      }
      return pn;
    }
    // Check if the given cell is outside of the map in i-the dimension
    inline bool isOutsideXYZ(const Veci<Dim> &n, int i) const
    {
      return n(i) < 0 || n(i) >= dim_(i);
    }
    // Check if the cell is free by index
    inline bool isFree(int idx) const
    {
      return map_[idx] == val_free_;
    }
    // Check if the cell is unknown by index
    inline bool isUnknown(int idx) const
    {
      return map_[idx] == val_unknown_;
    }
    // Check if the cell is occupied by index
    inline bool isOccupied(int idx) const
    {
      return map_[idx] > val_free_;
    }

    inline void setOccupied(const Veci<Dim> &pn)
    {
      int index = getIndex(pn);
      if (index >= 0 && index < total_size_)
      { // check that the point is inside the map
        map_[getIndex(pn)] = val_occ_;
      }
    }

    inline void setFree(const Veci<Dim> &pn)
    {
      int index = getIndex(pn);
      if (index >= 0 && index < total_size_)
      { // check that the point is inside the map
        map_[index] = val_free_;
      }
    }

    // set Free all the voxels that are in a 3d cube centered at center and with side/2=d
    inline void setFreeVoxelAndSurroundings(const Veci<Dim> &center, const float d)
    {
      int n_voxels = std::round(d / res_ + 0.5); // convert distance to number of voxels
      for (int ix = -n_voxels; ix <= n_voxels; ix++)
      {
        for (int iy = -n_voxels; iy <= n_voxels; iy++)
        {
          for (int iz = -n_voxels; iz <= n_voxels; iz++)
          {
            Veci<Dim> voxel = center + Veci<Dim>(ix, iy, iz); // Int coordinates of the voxel I'm going to clear

            // std::cout << "Clearing" << voxel.transpose() << std::endl;
            setFree(voxel);
          }
        }
      }
    }

    // Check if the cell is outside by coordinate
    inline bool isOutsideOldMap(const Veci<Dim> &pn) const
    {
      for (int i = 0; i < Dim; i++)
        if (pn(i) < 0 || pn(i) >= prev_dim_(i))
          return true;
      return false;
    }
    // Check if the cell is outside by coordinate
    inline bool isOutside(const Veci<Dim> &pn) const
    {
      for (int i = 0; i < Dim; i++)
        if (pn(i) < 0 || pn(i) >= dim_(i))
          return true;
      return false;
    }
    inline bool isOutside(const Vecf<Dim> &pt) const
    {
      return isOutside(floatToInt(pt));
    }
    // Check if the given cell is free by coordinate
    inline bool isFree(const Veci<Dim> &pn) const
    {
      if (isOutside(pn))
        return false;
      else
        return isFree(getIndex(pn));
    }
    inline bool isFree(const Vecf<Dim> &pt) const
    {
      return isFree(floatToInt(pt));
    }
    // Check if the given cell is occupied by coordinate
    inline bool isOccupied(const Veci<Dim> &pn) const
    {
      if (isOutside(pn))
        return false;
      else
        return isOccupied(getIndex(pn));
    }
    inline bool isOccupied(const Vecf<Dim> &pt) const
    {
      return isOccupied(floatToInt(pt));
    }
    inline bool isStaticOccupied(const Veci<Dim> &pn) const
    {
      if (isOutside(pn))
        return false;
      else
        return isStaticOccupied(getIndex(pn));
    }
    inline bool isStaticOccupied(const Vecf<Dim> &pt) const
    {
      return isStaticOccupied(floatToInt(pt));
    }
    // Check if the given cell is unknown by coordinate
    inline bool isUnknown(const Veci<Dim> &pn) const
    {
      if (isOutside(pn))
        return false;
      return map_[getIndex(pn)] == val_unknown_;
    }
    inline bool isUnknown(const Vecf<Dim> &pt) const
    {
      return isUnknown(floatToInt(pt));
    }

    // Print basic information about the util
    void info()
    {
      Vecf<Dim> range = dim_.template cast<decimal_t>() * res_;
      std::cout << "MapUtil Info ========================== " << std::endl;
      std::cout << "   res: [" << res_ << "]" << std::endl;
      std::cout << "   origin: [" << origin_d_.transpose() << "]" << std::endl;
      std::cout << "   range: [" << range.transpose() << "]" << std::endl;
      std::cout << "   dim: [" << dim_.transpose() << "]" << std::endl;
    };

    // Float position to discrete cell coordinate
    inline Veci<Dim> floatToInt(const Vecf<Dim> &pt) const
    {
      return ((pt - origin_d_) / res_ - Vecf<Dim>::Constant(0.5)).template cast<int>();
    }

    // Discrete cell coordinate to float position
    inline Vecf<Dim> intToFloat(const Veci<Dim> &pn) const
    {
      return (pn.template cast<decimal_t>() + Vecf<Dim>::Constant(0.5)) * res_ + origin_d_;
    }

    // Raytrace from float point pt1 to pt2
    inline vec_Veci<Dim> rayTrace(const Vecf<Dim> &pt1, const Vecf<Dim> &pt2) const
    {
      Vecf<Dim> diff = pt2 - pt1;
      decimal_t k = 0.1;
      int max_diff = (diff / res_).template lpNorm<Eigen::Infinity>() / k;
      decimal_t s = 1.0 / max_diff;
      Vecf<Dim> step = diff * s;

      vec_Veci<Dim> pns;
      Veci<Dim> prev_pn = Veci<Dim>::Constant(-1);
      for (int n = 1; n < max_diff; n++)
      {
        Vecf<Dim> pt = pt1 + step * n;
        Veci<Dim> new_pn = floatToInt(pt);
        if (isOutside(new_pn))
          break;
        if (new_pn != prev_pn)
          pns.push_back(new_pn);
        prev_pn = new_pn;
      }
      return pns;
    }

    // Check if the ray from p1 to p2 is occluded
    inline bool isBlocked(const Vecf<Dim> &p1, const Vecf<Dim> &p2, int8_t val = 100) const
    {
      vec_Veci<Dim> pns = rayTrace(p1, p2);
      for (const auto &pn : pns)
      {
        if (map_[getIndex(pn)] >= val)
          return true;
      }
      return false;
    }

    // Compute vicinities
    void computeVicinityMapInteger(const vec_Vecf<3> &path, const std::vector<float> &local_box_size, Veci<3> &min_point_int, Veci<3> &max_point_int) const
    {
      // 1. Compute the global bounding box around the path.
      Vecf<3> min_point_float = Vecf<3>::Constant(std::numeric_limits<float>::max());
      Vecf<3> max_point_float = Vecf<3>::Constant(std::numeric_limits<float>::lowest());

      for (const auto &point : path)
      {
        // Inflate the local box size (using factor 1.5 as in your example)
        Vecf<3> inflated_local_box_size(1.5 * local_box_size[0], 1.5 * local_box_size[1], 1.5 * local_box_size[2]);
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
    vec_Vecf<Dim> getCloud_(CheckFunc check, const Veci<3> &min_point_int, const Veci<3> &max_point_int) const
    {
      vec_Vecf<Dim> cloud;

      // Reserve an estimated size (optional, just to reduce reallocations)
      cloud.reserve(static_cast<size_t>((max_point_int - min_point_int).prod()));

      if (Dim == 3)
      {
        // Store thread-local clouds to avoid serialization in critical section
        std::vector<vec_Vecf<Dim>> thread_clouds;
        int num_threads = 1;

#pragma omp parallel
        {
#pragma omp single
          {
            num_threads = omp_get_num_threads();
            thread_clouds.resize(num_threads);
          }

          int tid = omp_get_thread_num();
          vec_Vecf<Dim> &local_cloud = thread_clouds[tid];

// Collapse the three nested loops into one for OpenMP
#pragma omp for collapse(3) nowait
          for (int i = min_point_int(0); i < max_point_int(0); ++i)
          {
            for (int j = min_point_int(1); j < max_point_int(1); ++j)
            {
              for (int k = min_point_int(2); k < max_point_int(2); ++k)
              {
                Veci<3> pti(i, j, k);
                // Use the provided check function
                if (check(pti) && !isOutside(pti))
                {
                  Vecf<3> ptf = intToFloat(pti);
                  local_cloud.push_back(ptf);
                }
              }
            }
          }
        }

        // Merge all thread-local clouds serially (no critical section bottleneck)
        size_t total_size = 0;
        for (const auto &tc : thread_clouds)
          total_size += tc.size();
        cloud.reserve(total_size);
        for (const auto &tc : thread_clouds)
          cloud.insert(cloud.end(), tc.begin(), tc.end());
      }
      else if (Dim == 2)
      {
        // Store thread-local clouds to avoid serialization in critical section
        std::vector<vec_Vecf<Dim>> thread_clouds;
        int num_threads = 1;

#pragma omp parallel
        {
#pragma omp single
          {
            num_threads = omp_get_num_threads();
            thread_clouds.resize(num_threads);
          }

          int tid = omp_get_thread_num();
          vec_Vecf<Dim> &local_cloud = thread_clouds[tid];

#pragma omp for collapse(2) nowait
          for (int i = min_point_int(0); i < max_point_int(0); ++i)
          {
            for (int j = min_point_int(1); j < max_point_int(1); ++j)
            {
              Veci<3> pti(i, j, 0);
              if (check(pti) && !isOutside(pti))
              {
                Vecf<3> ptf = intToFloat(pti);
                local_cloud.push_back(ptf);
              }
            }
          }
        }

        // Merge all thread-local clouds serially (no critical section bottleneck)
        size_t total_size = 0;
        for (const auto &tc : thread_clouds)
          total_size += tc.size();
        cloud.reserve(total_size);
        for (const auto &tc : thread_clouds)
          cloud.insert(cloud.end(), tc.begin(), tc.end());
      }

      return cloud;
    }

    // Cloud-getter
    template <typename CheckFunc>
    vec_Vecf<Dim> getCloud(CheckFunc check, const vec_Vecf<3> &path, const std::vector<float> local_box_size) const
    {
      vec_Vecf<Dim> cloud;

      // Compute vicinty map integer
      Veci<3> min_point_int, max_point_int;
      computeVicinityMapInteger(path, local_box_size, min_point_int, max_point_int);

      return getCloud_(check, min_point_int, max_point_int);
    }

    // Cloud-getter
    template <typename CheckFunc>
    vec_Vecf<Dim> getCloud(CheckFunc check) const
    {
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

    // Get occupied voxels (useful for convex decomposition)
    inline vec_Vecf<Dim> getOccupiedCloud(const vec_Vecf<3> &path, const std::vector<float> local_box_size) const
    {
      // Get cloud for occupied cells
      return getCloud([this](const Veci<Dim> &pti) -> bool
                      { return isOccupied(pti); }, path, local_box_size);
    }

    // Get occupied voxels for the entire map (useful for visualization)
    inline vec_Vecf<Dim> getOccupiedCloud() const
    {
      // Get cloud for occupied cells
      return getCloud([this](const Veci<Dim> &pti) -> bool
                      { return isOccupied(pti); });
    }

    // Get free voxels
    inline vec_Vecf<Dim> getFreeCloud(const vec_Vecf<3> &path, const std::vector<float> local_box_size) const
    {
      // Get cloud for free cells
      return getCloud([this](const Veci<Dim> &pti) -> bool
                      { return isFree(pti); }, path, local_box_size);
    }

    // Get free voxels for the entire map (useful for visualization)
    inline vec_Vecf<Dim> getFreeCloud() const
    {
      // Get cloud for free cells
      return getCloud([this](const Veci<Dim> &pti) -> bool
                      { return isFree(pti); });
    }

    // Get unknown voxels
    inline vec_Vecf<Dim> getUnknownCloud(const vec_Vecf<3> &path, const std::vector<float> local_box_size) const
    {
      // Get cloud for unknown cells
      return getCloud([this](const Veci<Dim> &pti) -> bool
                      { return isUnknown(pti); }, path, local_box_size);
    }

    // Get unknown voxels for the entire map (useful for visualization)
    inline vec_Vecf<Dim> getUnknownCloud() const
    {
      // Get cloud for unknown cells
      return getCloud([this](const Veci<Dim> &pti) -> bool
                      { return isUnknown(pti); });
    }

    // Get static occupied voxels with unknown as occupied for a specific region (useful for convex decomposition)
    inline vec_Vecf<Dim> getOccupiedCloudWithUnknownAsOccupied(const vec_Vecf<3> &path, const std::vector<float> local_box_size) const
    {
      // Get cloud for occupied cells
      return getCloud([this](const Veci<Dim> &pti) -> bool
                      { return isOccupied(pti) || isUnknown(pti); }, path, local_box_size);
    }

    // Get number of unknown cells
    int countUnknownCells() const
    {
      return std::count(map_.begin(), map_.end(), val_unknown_);
    }

    // Get number of occupied cells
    int countOccupiedCells() const
    {
      return std::count(map_.begin(), map_.end(), val_occ_);
    }

    // Get number of free cells
    int countFreeCells() const
    {
      return std::count(map_.begin(), map_.end(), val_free_);
    }

    // Get the total number of cells
    int getTotalNumCells() const
    {
      return total_size_;
    }

    // Map entity
    Tmap map_;

    // Heat map layer (soft costs)
    std::vector<float> heat_;
    bool dynamic_heat_enabled_{false};
    bool dynamic_as_occupied_current_{true};
    bool dynamic_as_occupied_future_{true};
    float heat_w_{0.0f};

    // Dynamic heat parameters
    float heat_alpha0_{1.0f};
    float heat_alpha1_{2.0f};
    int heat_p_{2};
    int heat_q_{2};
    float heat_tau_ratio_{0.5f};
    float heat_gamma_{0.0f};
    float heat_Hmax_{10.0f};
    float dyn_base_inflation_m_{0.5f};
    int heat_num_samples_{15};
    std::vector<vec_Vecf<3>> dyn_pred_samples_;
    std::vector<float> dyn_pred_times_;

    // Static heat parameters
    bool static_heat_enabled_{false};
    float static_heat_alpha_{2.0f};
    int static_heat_p_{2};
    float static_heat_Hmax_{50.0f};
    float static_heat_rmax_m_{1.0f};
    float static_heat_default_radius_m_{0.5f};
    bool static_heat_boundary_only_{true};
    bool static_heat_apply_on_unknown_{false};
    bool static_heat_exclude_dynamic_{true};
    std::function<float(const Eigen::Vector3f &)> static_heat_radius_fn_;

    // Static heat offset cache
    struct StaticHeatOff { int dx, dy, dz; float d_m; };
    mutable std::vector<StaticHeatOff> static_heat_off_;
    mutable int static_heat_off_Rcell_{-1};
    mutable float static_heat_off_res_{-1.0f};
    mutable float static_heat_off_rmax_m_{-1.0f};
    mutable std::mutex static_heat_mutex_;

    // Obstacle max velocity
    float obst_max_vel_{1.0f};

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

  private:
    void ensureStaticHeatOffsets(int Rcell) const;
    void computeDynamicHeat(const vec_Vecf<3> &obst_pos,
                           const vec_Vecf<3> &obst_bbox,
                           double traj_max_time);
    void computeStaticHeat();
  };

  // Template specialization for 3D MapUtil heat methods
  template <>
  inline vec_Vecf<3> MapUtil<3>::getHeatCloud(float threshold) const
  {
    vec_Vecf<3> cloud;
    if (heat_.empty())
      return cloud;

    for (int x = 0; x < dim_(0); ++x)
    {
      for (int y = 0; y < dim_(1); ++y)
      {
        for (int z = 0; z < dim_(2); ++z)
        {
          const int idx = x + y * dim_(0) + z * dim_xy_;
          if (idx >= 0 && idx < (int)heat_.size() && heat_[idx] > threshold)
          {
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

  template <>
  inline float MapUtil<3>::getMaxHeat() const
  {
    if (heat_.empty())
      return 0.0f;
    return *std::max_element(heat_.begin(), heat_.end());
  }

  // ensureStaticHeatOffsets implementation
  template <>
  inline void MapUtil<3>::ensureStaticHeatOffsets(int Rcell) const
  {
    std::lock_guard<std::mutex> lock(static_heat_mutex_);
    if (static_heat_off_Rcell_ == Rcell &&
        static_heat_off_res_ == (float)res_ &&
        std::fabs(static_heat_off_rmax_m_ - static_heat_rmax_m_) < 1e-6f)
    {
      return;
    }

    static_heat_off_Rcell_ = Rcell;
    static_heat_off_res_ = (float)res_;
    static_heat_off_rmax_m_ = static_heat_rmax_m_;

    static_heat_off_.clear();
    static_heat_off_.reserve((2 * Rcell + 1) * (2 * Rcell + 1) * (2 * Rcell + 1));

    for (int dx = -Rcell; dx <= Rcell; ++dx)
    {
      for (int dy = -Rcell; dy <= Rcell; ++dy)
      {
        for (int dz = -Rcell; dz <= Rcell; ++dz)
        {
          const float d_m = (float)res_ * std::sqrt(float(dx * dx + dy * dy + dz * dz));
          if (d_m > static_heat_rmax_m_)
            continue;
          static_heat_off_.push_back({dx, dy, dz, d_m});
        }
      }
    }
  }

  // computeDynamicHeat implementation
  template <>
  inline void MapUtil<3>::computeDynamicHeat(const vec_Vecf<3> &obst_pos,
                                              const vec_Vecf<3> &obst_bbox,
                                              double traj_max_time)
  {
    const float Th = std::max(0.0f, (float)traj_max_time);
    const float tau_w = std::max(1e-3f, heat_tau_ratio_ * std::max(1e-3f, Th));

    // Determine time samples
    std::vector<float> t_samples;
    if (!dyn_pred_times_.empty())
    {
      t_samples = dyn_pred_times_;
    }
    else
    {
      const int M = std::max(2, heat_num_samples_);
      t_samples.resize(M);
      for (int j = 0; j < M; ++j)
        t_samples[j] = (float)j * Th / (float)(M - 1);
    }

    const size_t K = obst_pos.size();
    if (K == 0)
      return;

    // Lightweight pow for small integer exponents
    auto pow_fast = [](float x, int p) -> float
    {
      x = std::max(0.0f, x);
      switch (p)
      {
      case 1:
        return x;
      case 2:
        return x * x;
      case 3:
        return x * x * x;
      case 4:
      {
        const float x2 = x * x;
        return x2 * x2;
      }
      default:
        return std::pow(x, (float)p);
      }
    };

    // Precompute obstacle centers, bbox half-extents, and reachable radii
    std::vector<Eigen::Vector3f> ck_list(K);
    std::vector<Eigen::Vector3f> hk_list(K);
    std::vector<float> Rreach_list(K);

    const float R0 = std::max(0.0f, dyn_base_inflation_m_);
    for (size_t k = 0; k < K; ++k)
    {
      ck_list[k] = obst_pos[k].cast<float>();

      // Get bbox half-extents for this obstacle
      float hx = 0.4f, hy = 0.4f, hz = 0.4f;
      if (k < obst_bbox.size())
      {
        hx = obst_bbox[k].x();
        hy = obst_bbox[k].y();
        hz = obst_bbox[k].z();
      }
      hk_list[k] = Eigen::Vector3f(hx, hy, hz);

      // Reachable radius: bbox extent + motion
      const float max_extent = std::max({hx, hy, hz});
      Rreach_list[k] = max_extent + (float)obst_max_vel_ * Th;
    }

    // Precompute per-time-sample tube radii and time-decay weights
    const size_t J = t_samples.size();
    std::vector<float> Rj(J), Wj(J);
    for (size_t j = 0; j < J; ++j)
    {
      const float tj = std::max(0.0f, t_samples[j]);
      Rj[j] = R0 + heat_gamma_ * tj;
      Wj[j] = std::exp(-tj / tau_w);
    }

    // Precompute predicted centers (fallback to ck if unavailable)
    std::vector<Eigen::Vector3f> cj_flat(K * J);
    for (size_t k = 0; k < K; ++k)
    {
      for (size_t j = 0; j < J; ++j)
      {
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
    for (int idx = 0; idx < total_size_; ++idx)
    {
      // Heat is only relevant for traversable cells
      if (map_[idx] > val_free_)
        continue;

      const int ix = idx % dim0;
      const int iy = (idx / dim0) % dim1;
      const int iz = idx / plane;

      const float xw = origin_d_(0) + (ix + 0.5f) * res_;
      const float yw = origin_d_(1) + (iy + 0.5f) * res_;
      const float zw = origin_d_(2) + (iz + 0.5f) * res_;

      float best = 0.0f;

      for (size_t k = 0; k < K; ++k)
      {
        // Base reachable radius (finite horizon)
        float Hbase = 0.0f;
        const float Rreach = Rreach_list[k];
        if (Rreach > 1e-6f)
        {
          const Eigen::Vector3f &ck = ck_list[k];
          const Eigen::Vector3f &hk = hk_list[k];

          // Compute distance from point to box
          const float dx_abs = std::abs(xw - ck.x());
          const float dy_abs = std::abs(yw - ck.y());
          const float dz_abs = std::abs(zw - ck.z());

          const float dx_box = std::max(0.0f, dx_abs - hk.x());
          const float dy_box = std::max(0.0f, dy_abs - hk.y());
          const float dz_box = std::max(0.0f, dz_abs - hk.z());

          const float d2 = dx_box * dx_box + dy_box * dy_box + dz_box * dz_box;
          const float R2 = Rreach * Rreach;

          if (d2 <= R2)
          {
            const float d = std::sqrt(std::max(0.0f, d2));
            const float u = std::min(1.0f, std::max(0.0f, d / Rreach));
            Hbase = heat_alpha0_ * pow_fast(1.0f - u, heat_p_);
          }
        }

        // Tube bonus (max over time)
        float tube_max = 0.0f;
        const Eigen::Vector3f *cj_ptr = &cj_flat[k * J];
        const Eigen::Vector3f &hk = hk_list[k];

        for (size_t j = 0; j < J; ++j)
        {
          const float R = Rj[j];
          if (R <= 1e-6f)
            continue;

          const Eigen::Vector3f &cj = cj_ptr[j];

          // Distance from point to box at predicted position
          const float dx_abs = std::abs(xw - cj.x());
          const float dy_abs = std::abs(yw - cj.y());
          const float dz_abs = std::abs(zw - cj.z());

          const float dx_box = std::max(0.0f, dx_abs - hk.x());
          const float dy_box = std::max(0.0f, dy_abs - hk.y());
          const float dz_box = std::max(0.0f, dz_abs - hk.z());

          const float d2 = dx_box * dx_box + dy_box * dy_box + dz_box * dz_box;
          const float R2 = R * R;

          if (d2 > R2)
            continue;

          const float d = std::sqrt(std::max(0.0f, d2));
          const float u = std::min(1.0f, std::max(0.0f, d / R));
          const float g = pow_fast(1.0f - u, heat_q_);
          tube_max = std::max(tube_max, Wj[j] * g);
        }

        float Hk = Hbase + heat_alpha1_ * tube_max;
        if (heat_Hmax_ > 0.0f)
          Hk = std::min(Hk, heat_Hmax_);

        best = std::max(best, Hk);
      }

      heat_[idx] = std::max(heat_[idx], best);
    }
  }

  // computeStaticHeat implementation
  template <>
  inline void MapUtil<3>::computeStaticHeat()
  {
    if (static_heat_alpha_ <= 0.0f || static_heat_rmax_m_ <= 1e-6f)
      return;

    // Precompute offsets up to rmax
    const int Rcell = int(std::ceil(static_heat_rmax_m_ / res_));
    ensureStaticHeatOffsets(Rcell);
    const auto &off = static_heat_off_;

    const int dimX = dim_(0);
    const int dimY = dim_(1);
    const int dimZ = dim_(2);

    auto idx3_local = [&](int x, int y, int z)
    {
      return size_t(x) + size_t(dimX) * size_t(y) + size_t(dimX) * size_t(dimY) * size_t(z);
    };

    // Collect seeds (boundary occupied voxels)
    std::vector<size_t> seeds;
    seeds.reserve(size_t(total_size_ / 50) + 1);

    const int nx[6] = {+1, -1, 0, 0, 0, 0};
    const int ny[6] = {0, 0, +1, -1, 0, 0};
    const int nz[6] = {0, 0, 0, 0, +1, -1};

    for (int z = 0; z < dimZ; ++z)
    {
      for (int y = 0; y < dimY; ++y)
      {
        for (int x = 0; x < dimX; ++x)
        {
          const size_t lin = idx3_local(x, y, z);
          if (map_[lin] != val_occ_)
            continue;

          if (!static_heat_boundary_only_)
          {
            seeds.push_back(lin);
            continue;
          }

          bool boundary = false;
          for (int k = 0; k < 6; ++k)
          {
            const int x2 = x + nx[k], y2 = y + ny[k], z2 = z + nz[k];
            if (x2 < 0 || x2 >= dimX || y2 < 0 || y2 >= dimY || z2 < 0 || z2 >= dimZ)
            {
              boundary = true;
              break;
            }
            const size_t nlin = idx3_local(x2, y2, z2);
            if (map_[nlin] != val_occ_)
            {
              boundary = true;
              break;
            }
          }

          if (boundary)
            seeds.push_back(lin);
        }
      }
    }

    // Apply halo (max aggregation) around each seed
    for (const auto &lin : seeds)
    {
      const int x0 = int(lin % size_t(dimX));
      const int y0 = int((lin / size_t(dimX)) % size_t(dimY));
      const int z0 = int(lin / (size_t(dimX) * size_t(dimY)));

      // Seed voxel center in world coordinates (for radius function)
      const float xc0 = origin_d_(0) + (x0 + 0.5f) * res_;
      const float yc0 = origin_d_(1) + (y0 + 0.5f) * res_;
      const float zc0 = origin_d_(2) + (z0 + 0.5f) * res_;
      const Eigen::Vector3f seed_world(xc0, yc0, zc0);

      float Rm = static_heat_default_radius_m_;
      if (static_heat_radius_fn_)
        Rm = static_heat_radius_fn_(seed_world);

      // Clamp radius for safety/perf
      Rm = std::clamp(Rm, 0.0f, static_heat_rmax_m_);
      if (Rm <= 1e-6f)
        continue;

      for (const auto &o : off)
      {
        if (o.d_m > Rm)
          continue;

        const int x = x0 + o.dx, y = y0 + o.dy, z = z0 + o.dz;
        if (x < 0 || x >= dimX || y < 0 || y >= dimY || z < 0 || z >= dimZ)
          continue;

        const size_t idx = idx3_local(x, y, z);

        // Never override hard obstacles
        if (map_[idx] > val_free_)
          continue;

        const float u = std::min(1.0f, std::max(0.0f, o.d_m / Rm));
        const float base = 1.0f - u;
        float power_result;
        if (static_heat_p_ == 2)
        {
          power_result = base * base;
        }
        else if (static_heat_p_ == 3)
        {
          power_result = base * base * base;
        }
        else if (static_heat_p_ == 4)
        {
          const float base2 = base * base;
          power_result = base2 * base2;
        }
        else
        {
          power_result = std::pow(base, float(static_heat_p_));
        }
        float w = static_heat_alpha_ * power_result;

        if (static_heat_Hmax_ > 0.0f)
          w = std::min(w, static_heat_Hmax_);

        if (w > 0.0f)
        {
          heat_[idx] = std::max(heat_[idx], w);
        }
      }
    }
  }

  typedef MapUtil<3> VoxelMapUtil;

} // namespace mighty

#endif
