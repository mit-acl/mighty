#pragma once

#include <cmath>
#include <memory>
#include <vector>

#include <Eigen/Dense>
#include <nav_msgs/msg/occupancy_grid.hpp>

/** @brief 2D Euclidean Signed Distance Field grid with O(1) bilinear queries.
 *
 *  Constructed from a nav_msgs::OccupancyGrid published by the ACL mapper.
 *  The mapper encodes linearly in distance: value = 100 * (1 - d / d_max),
 *  where d_max = truncation_distance * resolution (meters). This class
 *  pre-converts to float distances in meters for fast bilinear interpolation
 *  during optimization.
 *
 *  Instances are immutable (const after construction) and can be shared
 *  across threads via shared_ptr<const EsdfGrid2D>.
 */
class EsdfGrid2D {
 public:
  /** @brief Construct from an OccupancyGrid message.
   *  @param msg  The OccupancyGrid from the mapper's esdf_2d_topic.
   *  @param truncation_voxels  Truncation distance in voxels (must match mapper config).
   */
  static std::shared_ptr<const EsdfGrid2D> fromOccupancyGrid(
      const nav_msgs::msg::OccupancyGrid& msg, int truncation_voxels) {
    auto grid = std::make_shared<EsdfGrid2D>();
    grid->resolution_ = msg.info.resolution;
    // Mapper shifts origin by -0.5*resolution for rviz alignment; undo for queries
    grid->origin_x_ = msg.info.origin.position.x + 0.5 * grid->resolution_;
    grid->origin_y_ = msg.info.origin.position.y + 0.5 * grid->resolution_;
    grid->inv_resolution_ = 1.0 / grid->resolution_;
    grid->width_ = static_cast<int>(msg.info.width);
    grid->height_ = static_cast<int>(msg.info.height);

    grid->max_distance_ = static_cast<double>(truncation_voxels) * grid->resolution_;

    // Pre-compute float distance grid (meters) from int8 OccupancyGrid encoding.
    // Mapper encodes linearly: value = 100 * (1 - d / d_max), so
    //   d = d_max * (1 - value / 100).
    const size_t n = msg.data.size();
    grid->dist_.resize(n);
    for (size_t i = 0; i < n; ++i) {
      int8_t val = msg.data[i];
      if (val <= 0) {
        grid->dist_[i] = grid->max_distance_;  // far from obstacles
      } else if (val >= 100) {
        grid->dist_[i] = 0.0;  // on obstacle
      } else {
        grid->dist_[i] = grid->max_distance_ * (1.0 - static_cast<double>(val) / 100.0);
      }
    }
    return grid;
  }

  /** @brief Query interpolated distance at world coordinates (x, y).
   *  @return Distance in meters. Returns max_distance_ if out of bounds. */
  double queryDistance(double x, double y) const {
    // OccupancyGrid origin is corner of cell (0,0); cell centers are at origin + (i+0.5)*res.
    // Subtract 0.5 so that a query at cell center maps to integer grid index.
    double gx = (x - origin_x_) * inv_resolution_;
    double gy = (y - origin_y_) * inv_resolution_;
    int ix = static_cast<int>(std::floor(gx));
    int iy = static_cast<int>(std::floor(gy));

    // Out-of-bounds: return max distance (no penalty)
    if (ix < 0 || ix >= width_ - 1 || iy < 0 || iy >= height_ - 1)
      return max_distance_;

    double fx = gx - ix;
    double fy = gy - iy;

    double d00 = dist_[iy * width_ + ix];
    double d10 = dist_[iy * width_ + ix + 1];
    double d01 = dist_[(iy + 1) * width_ + ix];
    double d11 = dist_[(iy + 1) * width_ + ix + 1];

    return (1.0 - fx) * (1.0 - fy) * d00 + fx * (1.0 - fy) * d10 +
           (1.0 - fx) * fy * d01 + fx * fy * d11;
  }

  /** @brief Query analytical gradient [dd/dx, dd/dy] of the bilinear interpolation.
   *  @return 2D gradient in meters/meter. Returns (0,0) if out of bounds. */
  Eigen::Vector2d queryGradient(double x, double y) const {
    double gx = (x - origin_x_) * inv_resolution_;
    double gy = (y - origin_y_) * inv_resolution_;
    int ix = static_cast<int>(std::floor(gx));
    int iy = static_cast<int>(std::floor(gy));

    if (ix < 0 || ix >= width_ - 1 || iy < 0 || iy >= height_ - 1)
      return Eigen::Vector2d::Zero();

    double fx = gx - ix;
    double fy = gy - iy;

    double d00 = dist_[iy * width_ + ix];
    double d10 = dist_[iy * width_ + ix + 1];
    double d01 = dist_[(iy + 1) * width_ + ix];
    double d11 = dist_[(iy + 1) * width_ + ix + 1];

    // Analytical gradient of bilinear interpolation
    double dd_dx = inv_resolution_ * ((1.0 - fy) * (d10 - d00) + fy * (d11 - d01));
    double dd_dy = inv_resolution_ * ((1.0 - fx) * (d01 - d00) + fx * (d11 - d10));

    return Eigen::Vector2d(dd_dx, dd_dy);
  }

  /** @brief Check if world point is within the grid (with 1-cell margin for bilinear). */
  bool isInBounds(double x, double y) const {
    double gx = (x - origin_x_) * inv_resolution_;
    double gy = (y - origin_y_) * inv_resolution_;
    return gx >= 0.0 && gx < width_ - 1.0 && gy >= 0.0 && gy < height_ - 1.0;
  }

  double maxDistance() const { return max_distance_; }

 private:
  double origin_x_ = 0.0, origin_y_ = 0.0;
  double resolution_ = 0.2;
  double inv_resolution_ = 5.0;
  int width_ = 0, height_ = 0;
  double max_distance_ = 2.0;
  std::vector<double> dist_;  // pre-computed distance in meters, row-major
};
