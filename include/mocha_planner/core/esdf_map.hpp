#pragma once

/**
 * @file esdf_map.hpp
 * @brief 2D Euclidean Signed Distance Field using exact EDT
 *        (Felzenszwalb-Huttenlocher algorithm).
 *
 * No ROS dependencies. Input is raw grid data.
 * Provides sub-pixel distance and gradient queries via bilinear interpolation.
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <Eigen/Core>

namespace mocha {

class EsdfMap {
public:
  EsdfMap() = default;

  /**
   * @brief Build ESDF from a raw occupancy grid.
   * @param data        Row-major occupancy values (size = width * height)
   * @param width       Grid width in cells
   * @param height      Grid height in cells
   * @param resolution  Cell size in meters
   * @param origin_x    World X of grid origin
   * @param origin_y    World Y of grid origin
   * @param obstacle_threshold  Values >= this are obstacles
   * @param max_distance  Clamp distance to this value (meters)
   * @return true on success
   */
  bool build(const std::vector<int8_t>& data,
             int width, int height,
             double resolution,
             double origin_x, double origin_y,
             int8_t obstacle_threshold = 50,
             double max_distance = 5.0);

  /**
   * @brief Query distance and normalized gradient at a world position.
   * @param pos   World position (x, y)
   * @param dist  Output: signed distance (positive = free, negative = inside obstacle)
   * @param grad  Output: normalized gradient (points away from nearest obstacle)
   * @return true if query is valid
   */
  bool query(const Eigen::Vector2d& pos, double& dist, Eigen::Vector2d& grad) const;

  bool isReady() const { return ready_; }
  int width() const { return width_; }
  int height() const { return height_; }
  double resolution() const { return resolution_; }
  double originX() const { return origin_x_; }
  double originY() const { return origin_y_; }

  /// Direct access to raw SDF grid (for obstacle reference point extraction etc.)
  const std::vector<double>& sdfGrid() const { return sdf_; }

private:
  bool ready_{false};
  int width_{0}, height_{0};
  double resolution_{0.05};
  double origin_x_{0.0}, origin_y_{0.0};
  double max_distance_{5.0};

  std::vector<double> sdf_;
  std::vector<double> grad_x_;
  std::vector<double> grad_y_;

  /**
   * @brief 1D squared distance transform (Felzenszwalb-Huttenlocher).
   *
   * Computes D[q] = min_p { (q - p)^2 + f[p] } for q = 0..n-1
   * using the lower envelope of parabolas in O(n).
   *
   * @param f   Input array (0 for seed cells, INF for others)
   * @param n   Array length
   * @param d   Output: squared distance values
   */
  static void dt1d(const double* f, int n, double* d);

  /**
   * @brief 2D exact EDT via separable 1D passes.
   *
   * Input: binary grid (true = seed cell with distance 0).
   * Output: Euclidean distance for each cell.
   */
  static void edt2d(const std::vector<bool>& seeds, int w, int h,
                    double resolution, double max_dist,
                    std::vector<double>& dist_out);
};

} // namespace mocha
