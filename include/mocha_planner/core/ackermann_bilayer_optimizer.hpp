#pragma once

#include <vector>
#include <Eigen/Core>

#include "mocha_planner/core/types.hpp"

namespace mocha {

struct AckermannBilayerParameters
{
  double wheel_base{0.6};
  double v_max{2.0};
  double phi_max{0.55};
  double phi_dot_max{1.0};
  double a_lon_max{2.0};
  double a_lat_max{2.0};
  double min_pseudo_speed{0.3};
  double boundary_tangent_norm{1.0};

  double w_jerk{1.0};
  double w_time{10.0};
  double w_distance{50.0};
  double w_pseudo_speed{2000.0};
  double w_speed{1000.0};
  double w_steering{1000.0};
  double w_steering_rate{500.0};
  double w_longitudinal_acc{500.0};
  double w_lateral_acc{500.0};
  double w_obstacle{1500.0};
  double safe_distance{0.2};

  int kappa{12};
  int gear_sign{1};
  int n_segments{0};
  bool use_vmc{true};

  double start_yaw{0.0};
  double end_yaw{0.0};
  Eigen::Vector2d start_waypoint{Eigen::Vector2d::Zero()};
  Eigen::Vector2d end_waypoint{Eigen::Vector2d::Zero()};
  Eigen::Vector2d ref_point{Eigen::Vector2d::Zero()};

  std::vector<Eigen::Vector2d> prey_points;
  std::vector<Eigen::Vector2d> footprint_points;
  std::vector<double> initial_time_durations;
  std::vector<double> initial_pseudo_durations;
  EsdfQueryFunc esdf_query;
};

struct AckermannBilayerTrajectory
{
  Eigen::MatrixXd geom_coeffs;
  Eigen::MatrixXd time_coeffs;
  Eigen::VectorXd pseudo_durations;
  Eigen::VectorXd time_durations;
  std::vector<Eigen::Vector2d> waypoints;
  double total_duration{0.0};
  double total_pseudo_length{0.0};
  double final_cost{0.0};
  int gear_sign{1};

  bool isValid() const
  {
    return geom_coeffs.size() > 0 &&
           time_coeffs.size() > 0 &&
           pseudo_durations.size() > 0 &&
           time_durations.size() > 0;
  }
};

struct AckermannBodyState
{
  Eigen::Vector2d position{Eigen::Vector2d::Zero()};
  double yaw{0.0};
  double speed{0.0};
  double curvature{0.0};
  double steering{0.0};
  double steering_rate{0.0};
  double longitudinal_acc{0.0};
  double lateral_acc{0.0};
};

class AckermannBilayerOptimizer {
public:
  struct LbfgsData {
    const AckermannBilayerParameters* params{nullptr};
  };

  static bool optimize(const AckermannBilayerParameters& params,
                       AckermannBilayerTrajectory& trajectory);

  static double costFunctionVector(const Eigen::VectorXd& x_vec,
                                   Eigen::VectorXd* g_out,
                                   void* data);

  static AckermannBodyState evaluateBodyState(const AckermannBilayerTrajectory& trajectory,
                                              double wheel_base,
                                              double t);

  static std::vector<Eigen::Vector2d> samplePositions(const AckermannBilayerTrajectory& trajectory,
                                                      double dt);

  static Eigen::VectorXd forwardT(const Eigen::VectorXd& tau);
  static Eigen::VectorXd backwardT(const Eigen::VectorXd& T);
  static Eigen::VectorXd backwardGradT(const Eigen::VectorXd& tau,
                                       const Eigen::VectorXd& gradT);
};

}  // namespace mocha
