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
  using Tmap = std::vector<int>;
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
    MapUtil(float res, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, float inflation)
    {

      /* --------- Initialize parameters --------- */
      setInflation(inflation);                                                                                                                               // Set inflation
      setResolution(res);                                                                                                                                    // Set the resolution
      setMapSize(x_min, x_max, y_min, y_max, z_min, z_max);                                                                                                  // Set the cells and z_boundaries
    }

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
        double inflation)
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
      total_size_ = total;
      origin_d_ = origin;
      center_map_ = center_map;
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

        // Get the neighboring indices
        std::vector<int> neighbor_indices;

        // Increase the radius until a free point is found
        for (float radius = 1.0; radius < 5.0; radius += 0.5) // TODO: expose the radius as a parameter
        {
          neighbor_indices.clear();
          getNeighborIndices(point_int, neighbor_indices, radius);

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
      return Dim == 2 ? pn(0) + dim_(0) * pn(1) : pn(0) + dim_(0) * pn(1) + dim_(0) * dim_(1) * pn(2);
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
#pragma omp parallel
        {
          vec_Vecf<Dim> local_cloud;
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
#pragma omp critical
          {
            cloud.insert(cloud.end(), local_cloud.begin(), local_cloud.end());
          }
        }
      }
      else if (Dim == 2)
      {
#pragma omp parallel
        {
          vec_Vecf<Dim> local_cloud;
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
#pragma omp critical
          {
            cloud.insert(cloud.end(), local_cloud.begin(), local_cloud.end());
          }
        }
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
  };

  typedef MapUtil<3> VoxelMapUtil;

} // namespace mighty

#endif
