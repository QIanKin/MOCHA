#include "mocha_planner/core/trajectory_utils.hpp"
#include "mocha_planner/core/mco_optimizer.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>

namespace mocha {

namespace {

struct SegmentSample
{
  int index{0};
  double local_time{0.0};
};

SegmentSample locateSegmentSample(const Eigen::VectorXd& T, double t)
{
  const int n_segments = static_cast<int>(T.size());
  SegmentSample sample;
  if (n_segments <= 0) {
    return sample;
  }

  const double total_time = T.sum();
  const double clamped_t = std::max(0.0, std::min(t, total_time));

  double time_accum = 0.0;
  for (int seg = 0; seg < n_segments; ++seg) {
    const double seg_end = time_accum + T(seg);
    if (clamped_t < seg_end || seg == n_segments - 1) {
      sample.index = seg;
      sample.local_time = std::max(0.0, std::min(clamped_t - time_accum, T(seg)));
      return sample;
    }
    time_accum = seg_end;
  }

  sample.index = n_segments - 1;
  sample.local_time = T(n_segments - 1);
  return sample;
}

std::vector<Eigen::Vector2d> deduplicateWaypoints(const std::vector<Eigen::Vector2d>& input)
{
  std::vector<Eigen::Vector2d> output;
  output.reserve(input.size());
  for (const auto& point : input) {
    if (output.empty() || (point - output.back()).norm() > 1e-4) {
      output.push_back(point);
    }
  }
  return output;
}

bool solveReferenceTrajectory(
    const std::vector<Eigen::Vector2d>& waypoints,
    const Eigen::VectorXd& segment_times,
    const Eigen::Vector2d& start_vel,
    const Eigen::Vector2d& start_acc,
    const Eigen::Vector2d& end_vel,
    const Eigen::Vector2d& end_acc,
    McoTrajectory& trajectory)
{
  if (waypoints.size() < 2 || segment_times.size() != static_cast<int>(waypoints.size()) - 1) {
    return false;
  }

  McoParameters params;
  params.n_segments = static_cast<int>(waypoints.size()) - 1;
  params.start_waypoint = waypoints.front();
  params.end_waypoint = waypoints.back();
  params.start_vel = start_vel;
  params.start_acc = start_acc;
  params.end_vel = end_vel;
  params.end_acc = end_acc;

  Eigen::MatrixXd M_banded;
  Eigen::MatrixXd b;
  McoOptimizer::buildMcoMatrix(waypoints, segment_times, params, M_banded, b);
  if (M_banded.size() == 0 || b.size() == 0) {
    return false;
  }

  Eigen::MatrixXd factored = M_banded;
  Eigen::MatrixXd coeffs =
      McoOptimizer::solveBandedSystem(factored, b, params.n_coeffs, params.n_coeffs);
  if (!coeffs.allFinite()) {
    return false;
  }

  trajectory.coeffs = coeffs;
  trajectory.T = segment_times;
  trajectory.total_duration = segment_times.sum();
  trajectory.final_cost = 0.0;
  return true;
}

}  // namespace

// Evaluate position on the piecewise polynomial trajectory.
Eigen::Vector2d TrajectoryUtils::evalPosition(const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
                                               int n_coeffs, double t) {
  const SegmentSample sample = locateSegmentSample(T, t);
  int row0 = sample.index * n_coeffs;
  double px = 0.0, py = 0.0, tp = 1.0;
  for (int i = 0; i < n_coeffs; ++i) {
    px += coeffs(row0 + i, 0) * tp;
    py += coeffs(row0 + i, 1) * tp;
    tp *= sample.local_time;
  }
  return Eigen::Vector2d(px, py);
}

// Evaluate first derivative on the piecewise polynomial trajectory.
Eigen::Vector2d TrajectoryUtils::evalVelocity(const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
                                               int n_coeffs, double t) {
  const SegmentSample sample = locateSegmentSample(T, t);
  int row0 = sample.index * n_coeffs;
  double vx = 0.0, vy = 0.0, tp = 1.0;
  for (int i = 1; i < n_coeffs; ++i) {
    vx += i * coeffs(row0 + i, 0) * tp;
    vy += i * coeffs(row0 + i, 1) * tp;
    tp *= sample.local_time;
  }
  return Eigen::Vector2d(vx, vy);
}

// Evaluate second derivative on the piecewise polynomial trajectory.
Eigen::Vector2d TrajectoryUtils::evalAcceleration(const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
                                                   int n_coeffs, double t) {
  const SegmentSample sample = locateSegmentSample(T, t);
  int row0 = sample.index * n_coeffs;
  double ax = 0.0, ay = 0.0, tp = 1.0;
  for (int i = 2; i < n_coeffs; ++i) {
    ax += i * (i - 1) * coeffs(row0 + i, 0) * tp;
    ay += i * (i - 1) * coeffs(row0 + i, 1) * tp;
    tp *= sample.local_time;
  }
  return Eigen::Vector2d(ax, ay);
}

std::vector<Eigen::Vector2d> TrajectoryUtils::resamplePolyline(
    const std::vector<Eigen::Vector2d>& polyline, double step) {
  if (polyline.size() < 2 || step <= 0) return polyline;

  // Compute total length
  double total_length = 0.0;
  for (size_t i = 1; i < polyline.size(); ++i) {
    total_length += (polyline[i] - polyline[i-1]).norm();
  }
  if (total_length < 1e-9) return polyline;

  std::vector<Eigen::Vector2d> result;
  result.push_back(polyline.front());

  double next_sample = step;
  size_t seg = 0;
  double seg_start = 0.0;

  while (seg < polyline.size() - 1) {
    double seg_len = (polyline[seg+1] - polyline[seg]).norm();
    double seg_end = seg_start + seg_len;

    while (next_sample <= seg_end && next_sample < total_length - 1e-9) {
      double t = (next_sample - seg_start) / seg_len;
      Eigen::Vector2d p = polyline[seg] + t * (polyline[seg+1] - polyline[seg]);
      result.push_back(p);
      next_sample += step;
    }

    seg_start = seg_end;
    ++seg;
  }

  // Always include last point
  if ((result.back() - polyline.back()).norm() > 1e-9) {
    result.push_back(polyline.back());
  }

  return result;
}

std::vector<Eigen::Vector2d> TrajectoryUtils::sampleTrajectory(
    const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
    int n_coeffs, double dt) {
  std::vector<Eigen::Vector2d> points;
  double total = T.sum();
  for (double t = 0.0; t <= total; t += dt) {
    points.push_back(evalPosition(coeffs, T, n_coeffs, t));
  }
  // Ensure last point
  if (points.empty() || (points.back() - evalPosition(coeffs, T, n_coeffs, total)).norm() > 1e-9) {
    points.push_back(evalPosition(coeffs, T, n_coeffs, total));
  }
  return points;
}

double TrajectoryUtils::estimateMaxSpeed(
    const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
    int n_coeffs, int samples_per_segment) {
  if (T.size() == 0 || coeffs.size() == 0) {
    return 0.0;
  }

  const int clamped_samples = std::max(2, samples_per_segment);
  double peak_speed = 0.0;
  double time_accum = 0.0;
  for (int seg = 0; seg < T.size(); ++seg) {
    const double seg_duration = T(seg);
    for (int i = 0; i <= clamped_samples; ++i) {
      const double alpha = static_cast<double>(i) / static_cast<double>(clamped_samples);
      const double t = std::min(T.sum(), time_accum + alpha * seg_duration);
      peak_speed = std::max(peak_speed, evalVelocity(coeffs, T, n_coeffs, t).norm());
    }
    time_accum += seg_duration;
  }
  return peak_speed;
}

bool TrajectoryUtils::fitMinJerkTrajectory(
    const std::vector<Eigen::Vector2d>& waypoints,
    const Eigen::Vector2d& start_vel,
    const Eigen::Vector2d& start_acc,
    const Eigen::Vector2d& end_vel,
    const Eigen::Vector2d& end_acc,
    double max_speed,
    McoTrajectory& trajectory) {
  const std::vector<Eigen::Vector2d> unique_waypoints = deduplicateWaypoints(waypoints);
  if (unique_waypoints.size() < 2) {
    return false;
  }

  double desired_speed = std::max(0.2, max_speed / 1.5);
  for (int attempt = 0; attempt < 4; ++attempt) {
    Eigen::VectorXd segment_times(static_cast<int>(unique_waypoints.size()) - 1);
    for (int i = 0; i < segment_times.size(); ++i) {
      const double seg_length = (unique_waypoints[static_cast<size_t>(i) + 1] -
                                 unique_waypoints[static_cast<size_t>(i)]).norm();
      segment_times(i) = std::max(0.2, seg_length / std::max(0.1, desired_speed));
    }

    if (!solveReferenceTrajectory(
            unique_waypoints, segment_times,
            start_vel, start_acc, end_vel, end_acc,
            trajectory)) {
      return false;
    }

    const double peak_speed = estimateMaxSpeed(
        trajectory.coeffs, trajectory.T, 6, 12);
    if (peak_speed <= std::max(0.1, max_speed * 1.02) || desired_speed <= 0.15) {
      return true;
    }

    desired_speed /= 1.5;
  }

  return trajectory.isValid();
}

bool TrajectoryUtils::checkTrajectoryCollision(
    const Eigen::MatrixXd& coeffs, const Eigen::VectorXd& T,
    int n_coeffs, double dt,
    const EsdfQueryFunc& esdf_query,
    double safe_distance) {

  if (!esdf_query || T.size() == 0 || coeffs.size() == 0) {
    return false;
  }

  const double total_duration = T.sum();
  const int num_samples = std::max(10, static_cast<int>(std::ceil(total_duration / dt)));

  for (int i = 0; i <= num_samples; ++i) {
    const double t = total_duration * static_cast<double>(i) / static_cast<double>(num_samples);
    const Eigen::Vector2d pos = evalPosition(coeffs, T, n_coeffs, t);

    double dist = 0.0;
    Eigen::Vector2d grad;
    if (!esdf_query(pos, dist, grad)) {
      // Query failed, assume collision
      return true;
    }

    if (dist < safe_distance) {
      return true;
    }
  }

  return false;
}

} // namespace mocha
