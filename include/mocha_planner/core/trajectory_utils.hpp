#pragma once
#include <vector>
#include <Eigen/Dense>
#include "mocha_planner/core/types.hpp"

namespace mocha {

class TrajectoryUtils {
public:
  // Evaluate polynomial trajectory at time t
  // coeffs: (n_coeffs*n_segments) x dims matrix
  // T: segment durations vector
  static Eigen::Vector2d evalPosition(const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
                                       int n_coeffs, double t);
  static Eigen::Vector2d evalVelocity(const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
                                       int n_coeffs, double t);
  static Eigen::Vector2d evalAcceleration(const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
                                           int n_coeffs, double t);

  // Resample polyline to uniform spacing
  static std::vector<Eigen::Vector2d> resamplePolyline(
      const std::vector<Eigen::Vector2d>& polyline, double step);

  // Sample trajectory to discrete points at given time interval
  static std::vector<Eigen::Vector2d> sampleTrajectory(
      const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
      int n_coeffs, double dt);

  // Solve a continuous quintic reference trajectory through the given waypoints.
  // This routine only allocates segment times heuristically and enforces endpoint PVA.
  static bool fitMinJerkTrajectory(
      const std::vector<Eigen::Vector2d>& waypoints,
      const Eigen::Vector2d& start_vel,
      const Eigen::Vector2d& start_acc,
      const Eigen::Vector2d& end_vel,
      const Eigen::Vector2d& end_acc,
      double max_speed,
      McoTrajectory& trajectory);

  // Estimate the peak speed by dense sampling.
  static double estimateMaxSpeed(
      const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
      int n_coeffs, int samples_per_segment);

  // Check if a trajectory collides with obstacles using ESDF
  static bool checkTrajectoryCollision(
      const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
      int n_coeffs, double dt,
      const EsdfQueryFunc& esdf_query,
      double safe_distance);
};

} // namespace mocha
