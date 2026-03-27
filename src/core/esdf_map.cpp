#include "mocha_planner/core/esdf_map.hpp"

#include <queue>
#include <cstring>

namespace mocha {

// ============================================================
// Felzenszwalb-Huttenlocher 1D squared distance transform
// D[q] = min_p { (q - p)^2 + f[p] }
// Uses lower envelope of parabolas, O(n).
// ============================================================
void EsdfMap::dt1d(const double* f, int n, double* d)
{
  std::vector<int>    v(n);     // parabola vertex locations
  std::vector<double> z(n + 1); // intersection points between consecutive parabolas

  int k = 0;       // index of rightmost parabola in lower envelope
  v[0] = 0;
  z[0] = -1e20;
  z[1] =  1e20;

  for (int q = 1; q < n; ++q) {
    // intersection of parabola at v[k] and parabola at q
    double s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k]))
               / (2.0 * q - 2.0 * v[k]);

    while (s <= z[k]) {
      --k;
      s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k]))
          / (2.0 * q - 2.0 * v[k]);
    }
    ++k;
    v[k] = q;
    z[k] = s;
    z[k + 1] = 1e20;
  }

  k = 0;
  for (int q = 0; q < n; ++q) {
    while (z[k + 1] < (double)q) ++k;
    double dq = (double)(q - v[k]);
    d[q] = dq * dq + f[v[k]];
  }
}

// ============================================================
// 2D exact EDT via two separable 1D passes
// ============================================================
void EsdfMap::edt2d(const std::vector<bool>& seeds, int w, int h,
                    double resolution, double max_dist,
                    std::vector<double>& dist_out)
{
  const double INF = 1e20;
  const int N = w * h;

  // Intermediate buffer: squared distances in grid units
  std::vector<double> grid_sq(N);

  // Initialize: 0 for seed cells, INF for others
  for (int i = 0; i < N; ++i) {
    grid_sq[i] = seeds[i] ? 0.0 : INF;
  }

  // Pass 1: transform each row (horizontal)
  {
    std::vector<double> f(std::max(w, h));
    std::vector<double> d(std::max(w, h));
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) f[x] = grid_sq[y * w + x];
      dt1d(f.data(), w, d.data());
      for (int x = 0; x < w; ++x) grid_sq[y * w + x] = d[x];
    }

    // Pass 2: transform each column (vertical)
    for (int x = 0; x < w; ++x) {
      for (int y = 0; y < h; ++y) f[y] = grid_sq[y * w + x];
      dt1d(f.data(), h, d.data());
      for (int y = 0; y < h; ++y) grid_sq[y * w + x] = d[y];
    }
  }

  // Convert squared grid distances to metric distances, clamped
  double max_grid_dist = max_dist / resolution;
  double max_grid_sq = max_grid_dist * max_grid_dist;

  dist_out.resize(N);
  for (int i = 0; i < N; ++i) {
    double sq = std::min(grid_sq[i], max_grid_sq);
    dist_out[i] = std::sqrt(sq) * resolution;
  }
}

// ============================================================
// Build ESDF from raw occupancy grid
// ============================================================
bool EsdfMap::build(const std::vector<int8_t>& data,
                    int width, int height,
                    double resolution,
                    double origin_x, double origin_y,
                    int8_t obstacle_threshold,
                    double max_distance)
{
  width_ = width;
  height_ = height;
  resolution_ = resolution;
  origin_x_ = origin_x;
  origin_y_ = origin_y;
  max_distance_ = max_distance;

  if (width_ <= 0 || height_ <= 0 || resolution_ <= 0.0) return false;
  if ((int)data.size() < width_ * height_) return false;

  const int N = width_ * height_;

  // Classify cells
  std::vector<bool> is_obstacle(N, false);
  for (int i = 0; i < N; ++i) {
    is_obstacle[i] = (data[i] >= obstacle_threshold);
  }

  // EDT for free cells: seeds = obstacle cells, distance = how far free cells are from obstacles
  std::vector<double> dist_free;
  edt2d(is_obstacle, width_, height_, resolution_, max_distance_, dist_free);

  // EDT for obstacle cells: seeds = free cells, distance = how far obstacle cells are from free space
  std::vector<bool> is_free(N);
  for (int i = 0; i < N; ++i) is_free[i] = !is_obstacle[i];
  std::vector<double> dist_obst;
  edt2d(is_free, width_, height_, resolution_, max_distance_, dist_obst);

  // Combine into SDF: positive outside obstacles, negative inside
  sdf_.resize(N);
  for (int i = 0; i < N; ++i) {
    double d = is_obstacle[i] ? -dist_obst[i] : dist_free[i];
    sdf_[i] = std::max(-max_distance_, std::min(max_distance_, d));
  }

  // Compute gradients via central differences + normalization
  grad_x_.resize(N, 0.0);
  grad_y_.resize(N, 0.0);
  for (int y = 0; y < height_; ++y) {
    for (int x = 0; x < width_; ++x) {
      int idx = y * width_ + x;
      double gx = 0.0, gy = 0.0;

      // Central difference for x
      if (x > 0 && x < width_ - 1) {
        gx = (sdf_[idx + 1] - sdf_[idx - 1]) / (2.0 * resolution_);
      } else if (x == 0 && width_ > 1) {
        gx = (sdf_[idx + 1] - sdf_[idx]) / resolution_;
      } else if (x == width_ - 1 && width_ > 1) {
        gx = (sdf_[idx] - sdf_[idx - 1]) / resolution_;
      }

      // Central difference for y
      if (y > 0 && y < height_ - 1) {
        gy = (sdf_[(y + 1) * width_ + x] - sdf_[(y - 1) * width_ + x]) / (2.0 * resolution_);
      } else if (y == 0 && height_ > 1) {
        gy = (sdf_[(y + 1) * width_ + x] - sdf_[idx]) / resolution_;
      } else if (y == height_ - 1 && height_ > 1) {
        gy = (sdf_[idx] - sdf_[(y - 1) * width_ + x]) / resolution_;
      }

      // Normalize: exact ESDF satisfies ‖∇d‖ = 1 (Eikonal equation)
      double norm = std::sqrt(gx * gx + gy * gy);
      if (norm > 1e-6) {
        gx /= norm;
        gy /= norm;
      }

      grad_x_[idx] = gx;
      grad_y_[idx] = gy;
    }
  }

  ready_ = true;
  return true;
}

// ============================================================
// Query: bilinear interpolation of distance and gradient
// ============================================================
bool EsdfMap::query(const Eigen::Vector2d& pos, double& dist, Eigen::Vector2d& grad) const
{
  if (!ready_) return false;

  // Convert world position to grid coordinates (cell-center aligned)
  double u = (pos.x() - origin_x_) / resolution_ - 0.5;
  double v = (pos.y() - origin_y_) / resolution_ - 0.5;

  // Out of bounds: return max distance, zero gradient
  if (u < 0.0 || v < 0.0 || u >= width_ - 1.0 || v >= height_ - 1.0) {
    dist = max_distance_;
    grad.setZero();
    return true;
  }

  int ix = static_cast<int>(u);
  int iy = static_cast<int>(v);
  double fx = u - ix;
  double fy = v - iy;

  ix = std::max(0, std::min(ix, width_ - 2));
  iy = std::max(0, std::min(iy, height_ - 2));

  // Bilinear interpolation helper
  auto interp = [&](const std::vector<double>& field) -> double {
    double v00 = field[iy * width_ + ix];
    double v10 = field[iy * width_ + ix + 1];
    double v01 = field[(iy + 1) * width_ + ix];
    double v11 = field[(iy + 1) * width_ + ix + 1];
    return (1.0 - fx) * (1.0 - fy) * v00 + fx * (1.0 - fy) * v10 +
           (1.0 - fx) * fy * v01 + fx * fy * v11;
  };

  dist = interp(sdf_);
  grad.x() = interp(grad_x_);
  grad.y() = interp(grad_y_);

  // Re-normalize after interpolation (interpolation may slightly denormalize)
  double gn = grad.norm();
  if (gn > 1e-6) {
    grad /= gn;
  } else {
    grad.setZero();
  }

  return true;
}

} // namespace mocha
