#pragma once
#include <vector>
#include <string>
#include <Eigen/Dense>

namespace mocha {

struct ReferenceSelectionParams {
  double distance_ratio = 0.15;
  std::string manual_side = "auto";
};

class ReferenceSelector {
public:
  static Eigen::Vector2d select(
      const std::vector<Eigen::Vector2d>& prey_path,
      const ReferenceSelectionParams& params);
};

} // namespace mocha
