#include "mocha_planner/core/reference_selector.hpp"
#include <cmath>

namespace mocha {

Eigen::Vector2d ReferenceSelector::select(
    const std::vector<Eigen::Vector2d>& prey_path,
    const ReferenceSelectionParams& params) {

  if (prey_path.size() < 2) return Eigen::Vector2d::Zero();

  const Eigen::Vector2d& p_start = prey_path.front();
  const Eigen::Vector2d& p_end = prey_path.back();
  Eigen::Vector2d path_vec = p_end - p_start;
  double path_len = path_vec.norm();
  if (path_len < 1e-9) return Eigen::Vector2d::Zero();

  Eigen::Vector2d path_dir = path_vec / path_len;

  Eigen::Vector2d perp(-path_dir.y(), path_dir.x());
  bool select_left = true;
  if (params.manual_side == "left") {
    select_left = true;
  } else if (params.manual_side == "right") {
    select_left = false;
  } else {
    double signed_offset_sum = 0.0;
    double weight_sum = 0.0;
    for (size_t i = 1; i + 1 < prey_path.size(); ++i) {
      const double signed_offset = perp.dot(prey_path[i] - p_start);
      const double weight = 1.0;
      signed_offset_sum += weight * signed_offset;
      weight_sum += weight;
    }
    if (weight_sum > 1e-9) {
      select_left = (signed_offset_sum / weight_sum) >= 0.0;
    }
  }

  if (!select_left) {
    perp = -perp;
  }

  // Perpendicular offset
  Eigen::Vector2d mid = p_start + 0.5 * path_vec;
  double offset = params.distance_ratio * path_len;
  return mid + offset * perp;
}

} // namespace mocha
