#include "mocha_planner/core/ackermann_bilayer_optimizer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "mocha_planner/core/mco_optimizer.hpp"
#include "lbfgs.hpp"

namespace mocha {

namespace {

constexpr int kQuinticCoeffCount = 6;
constexpr int kQuinticSmoothOrder = 3;
constexpr double kHugeCost = 1.0e18;
constexpr double kMinSafeNorm = 1.0e-8;
constexpr double kDistanceSmoothEps = 1.0e-6;

struct SegmentSample
{
  int index{0};
  double local_time{0.0};
};

struct BilayerContext
{
  Eigen::VectorXd dT;
  Eigen::VectorXd dS;
  std::vector<Eigen::Vector2d> waypoints;

  Eigen::MatrixXd geom_factor;
  Eigen::MatrixXd geom_coeffs;

  Eigen::MatrixXd scale_factor;
  Eigen::MatrixXd scale_coeffs_dummy;
};

struct Poly2DEval
{
  Eigen::Vector2d p{Eigen::Vector2d::Zero()};
  Eigen::Vector2d d1{Eigen::Vector2d::Zero()};
  Eigen::Vector2d d2{Eigen::Vector2d::Zero()};
  Eigen::Vector2d d3{Eigen::Vector2d::Zero()};
  Eigen::Vector2d d4{Eigen::Vector2d::Zero()};
};

struct Poly1DEval
{
  double p{0.0};
  double d1{0.0};
  double d2{0.0};
  double d3{0.0};
  double d4{0.0};
};

struct MincoBackprop
{
  Eigen::MatrixXd grad_waypoints;
  Eigen::MatrixXd grad_head_pva;
  Eigen::MatrixXd grad_tail_pva;
  Eigen::VectorXd grad_durations;
};

Eigen::Vector2d headingVector(const double yaw)
{
  return Eigen::Vector2d(std::cos(yaw), std::sin(yaw));
}

void positiveCubicPenalty(const double violation, double& penalty, double& dpenalty)
{
  if (violation <= 0.0) {
    penalty = 0.0;
    dpenalty = 0.0;
    return;
  }
  penalty = violation * violation * violation;
  dpenalty = 3.0 * violation * violation;
}

std::vector<Eigen::Vector2d> buildVmcWaypoints(const Eigen::VectorXd& v_vec,
                                               const AckermannBilayerParameters& params)
{
  std::vector<Eigen::Vector2d> waypoints;
  waypoints.reserve(static_cast<size_t>(params.n_segments) + 1U);
  waypoints.push_back(params.start_waypoint);
  for (int i = 0; i < params.n_segments - 1; ++i) {
    const Eigen::Vector2d prey = params.prey_points[static_cast<size_t>(i) + 1U];
    waypoints.push_back(params.ref_point + v_vec(i) * (prey - params.ref_point));
  }
  waypoints.push_back(params.end_waypoint);
  return waypoints;
}

std::vector<Eigen::Vector2d> buildDirectWaypoints(const Eigen::VectorXd& q_vec,
                                                  const AckermannBilayerParameters& params)
{
  std::vector<Eigen::Vector2d> waypoints;
  waypoints.reserve(static_cast<size_t>(params.n_segments) + 1U);
  waypoints.push_back(params.start_waypoint);
  for (int i = 0; i < params.n_segments - 1; ++i) {
    waypoints.emplace_back(q_vec(2 * i), q_vec(2 * i + 1));
  }
  waypoints.push_back(params.end_waypoint);
  return waypoints;
}

void buildDistanceGradient(const std::vector<Eigen::Vector2d>& waypoints,
                           const double weight,
                           Eigen::MatrixXd& grad_q,
                           double& distance_cost)
{
  const int n_intermediate = std::max(0, static_cast<int>(waypoints.size()) - 2);
  grad_q.setZero(n_intermediate, 2);
  distance_cost = 0.0;

  if (weight <= 0.0 || waypoints.size() < 2) {
    return;
  }

  for (int i = 0; i + 1 < static_cast<int>(waypoints.size()); ++i) {
    const Eigen::Vector2d delta = waypoints[static_cast<size_t>(i) + 1U] - waypoints[static_cast<size_t>(i)];
    const double dist = std::sqrt(delta.squaredNorm() + kDistanceSmoothEps);
    const Eigen::Vector2d grad = delta / dist;
    distance_cost += dist;
    if (i > 0) {
      grad_q.row(i - 1) -= grad.transpose();
    }
    if (i < n_intermediate) {
      grad_q.row(i) += grad.transpose();
    }
  }

  distance_cost *= weight;
  grad_q *= weight;
}

Eigen::VectorXd backpropagateGradQtoSpatial(const Eigen::MatrixXd& grad_q,
                                            const AckermannBilayerParameters& params)
{
  const int n_intermediate = std::max(0, params.n_segments - 1);
  if (params.use_vmc) {
    Eigen::VectorXd grad_v = Eigen::VectorXd::Zero(n_intermediate);
    for (int i = 0; i < n_intermediate; ++i) {
      const Eigen::Vector2d prey = params.prey_points[static_cast<size_t>(i) + 1U];
      const Eigen::Vector2d dq_dv = prey - params.ref_point;
      grad_v(i) = grad_q.row(i).dot(dq_dv);
    }
    return grad_v;
  }

  Eigen::VectorXd grad_direct = Eigen::VectorXd::Zero(2 * n_intermediate);
  for (int i = 0; i < n_intermediate; ++i) {
    grad_direct(2 * i) = grad_q(i, 0);
    grad_direct(2 * i + 1) = grad_q(i, 1);
  }
  return grad_direct;
}

Eigen::Vector2d eval2DValue(const Eigen::MatrixXd& coeffs,
                            const int seg,
                            const double local_value,
                            const int derivative)
{
  const Eigen::VectorXd basis = McoOptimizer::getPolyBasis(local_value, 5, derivative);
  return coeffs.block(seg * kQuinticCoeffCount, 0, kQuinticCoeffCount, 2).transpose() * basis;
}

double eval1DValue(const Eigen::MatrixXd& coeffs_dummy,
                   const int seg,
                   const double local_value,
                   const int derivative)
{
  const Eigen::VectorXd basis = McoOptimizer::getPolyBasis(local_value, 5, derivative);
  return coeffs_dummy.block(seg * kQuinticCoeffCount, 0, kQuinticCoeffCount, 1).col(0).dot(basis);
}

Poly2DEval evalGeometrySegment(const Eigen::MatrixXd& coeffs,
                               const int seg,
                               const double local_s)
{
  Poly2DEval out;
  out.p = eval2DValue(coeffs, seg, local_s, 0);
  out.d1 = eval2DValue(coeffs, seg, local_s, 1);
  out.d2 = eval2DValue(coeffs, seg, local_s, 2);
  out.d3 = eval2DValue(coeffs, seg, local_s, 3);
  out.d4 = eval2DValue(coeffs, seg, local_s, 4);
  return out;
}

Poly1DEval evalScaleSegment(const Eigen::MatrixXd& coeffs_dummy,
                            const int seg,
                            const double local_t)
{
  Poly1DEval out;
  out.p = eval1DValue(coeffs_dummy, seg, local_t, 0);
  out.d1 = eval1DValue(coeffs_dummy, seg, local_t, 1);
  out.d2 = eval1DValue(coeffs_dummy, seg, local_t, 2);
  out.d3 = eval1DValue(coeffs_dummy, seg, local_t, 3);
  out.d4 = eval1DValue(coeffs_dummy, seg, local_t, 4);
  return out;
}

SegmentSample locateSegmentSample(const Eigen::VectorXd& durations, const double t)
{
  SegmentSample sample;
  const int n = static_cast<int>(durations.size());
  if (n <= 0) {
    return sample;
  }

  const double total = durations.sum();
  const double clamped = std::max(0.0, std::min(t, total));
  double acc = 0.0;
  for (int i = 0; i < n; ++i) {
    const double end = acc + durations(i);
    if (clamped < end || i == n - 1) {
      sample.index = i;
      sample.local_time = std::max(0.0, std::min(clamped - acc, durations(i)));
      return sample;
    }
    acc = end;
  }

  sample.index = n - 1;
  sample.local_time = durations(n - 1);
  return sample;
}

bool solveBilayerTrajectory(const AckermannBilayerParameters& params,
                            const Eigen::VectorXd& x_vec,
                            BilayerContext& context)
{
  const int n_segments = params.n_segments;
  const int n_intermediate = std::max(0, n_segments - 1);
  const int spatial_len = params.use_vmc ? n_intermediate : 2 * n_intermediate;

  if (x_vec.size() != 2 * n_segments + spatial_len) {
    return false;
  }

  const Eigen::VectorXd tauT = x_vec.head(n_segments);
  const Eigen::VectorXd tauS = x_vec.segment(n_segments, n_segments);
  const Eigen::VectorXd spatial =
      spatial_len > 0 ? x_vec.tail(spatial_len) : Eigen::VectorXd();

  context.dT = AckermannBilayerOptimizer::forwardT(tauT);
  context.dS = AckermannBilayerOptimizer::forwardT(tauS);
  if (!context.dT.allFinite() || !context.dS.allFinite()) {
    return false;
  }

  context.waypoints = spatial_len > 0
      ? (params.use_vmc ? buildVmcWaypoints(spatial, params)
                        : buildDirectWaypoints(spatial, params))
      : std::vector<Eigen::Vector2d>{params.start_waypoint, params.end_waypoint};

  if (static_cast<int>(context.waypoints.size()) != n_segments + 1) {
    return false;
  }

  const double q_boundary = params.boundary_tangent_norm;
  if (q_boundary <= params.min_pseudo_speed || params.min_pseudo_speed <= 0.0) {
    return false;
  }

  McoParameters geom_params;
  geom_params.n_segments = n_segments;
  geom_params.start_waypoint = params.start_waypoint;
  geom_params.end_waypoint = params.end_waypoint;
  geom_params.start_vel = params.gear_sign * q_boundary * headingVector(params.start_yaw);
  geom_params.end_vel = params.gear_sign * q_boundary * headingVector(params.end_yaw);
  geom_params.start_acc = Eigen::Vector2d::Zero();
  geom_params.end_acc = Eigen::Vector2d::Zero();

  Eigen::MatrixXd geom_b;
  McoOptimizer::buildMcoMatrix(context.waypoints, context.dS, geom_params, context.geom_factor, geom_b);
  if (context.geom_factor.size() == 0 || geom_b.size() == 0) {
    return false;
  }
  context.geom_coeffs = McoOptimizer::solveBandedSystem(
      context.geom_factor, geom_b, kQuinticCoeffCount, kQuinticCoeffCount);
  if (!context.geom_coeffs.allFinite()) {
    return false;
  }

  std::vector<Eigen::Vector2d> scale_waypoints;
  scale_waypoints.reserve(static_cast<size_t>(n_segments) + 1U);
  scale_waypoints.emplace_back(0.0, 0.0);
  double cumulative_s = 0.0;
  for (int i = 0; i < n_segments - 1; ++i) {
    cumulative_s += context.dS(i);
    scale_waypoints.emplace_back(cumulative_s, 0.0);
  }
  cumulative_s += context.dS(n_segments - 1);
  scale_waypoints.emplace_back(cumulative_s, 0.0);

  McoParameters scale_params;
  scale_params.n_segments = n_segments;
  scale_params.start_waypoint = Eigen::Vector2d::Zero();
  scale_params.end_waypoint = Eigen::Vector2d(cumulative_s, 0.0);
  scale_params.start_vel = Eigen::Vector2d::Zero();
  scale_params.end_vel = Eigen::Vector2d::Zero();
  scale_params.start_acc = Eigen::Vector2d::Zero();
  scale_params.end_acc = Eigen::Vector2d::Zero();

  Eigen::MatrixXd scale_b;
  McoOptimizer::buildMcoMatrix(scale_waypoints, context.dT, scale_params, context.scale_factor, scale_b);
  if (context.scale_factor.size() == 0 || scale_b.size() == 0) {
    return false;
  }
  context.scale_coeffs_dummy = McoOptimizer::solveBandedSystem(
      context.scale_factor, scale_b, kQuinticCoeffCount, kQuinticCoeffCount);
  if (!context.scale_coeffs_dummy.allFinite()) {
    return false;
  }

  return true;
}

MincoBackprop backpropagateMinco(const Eigen::VectorXd& durations,
                                 const Eigen::MatrixXd& total_grad_c,
                                 const Eigen::VectorXd& total_grad_t_direct,
                                 const Eigen::MatrixXd& factorized_m,
                                 const Eigen::MatrixXd& coeffs)
{
  MincoBackprop result;
  const int n_segments = static_cast<int>(durations.size());
  const int dims = static_cast<int>(total_grad_c.cols());
  const int s = kQuinticSmoothOrder;

  result.grad_waypoints.setZero(std::max(0, n_segments - 1), dims);
  result.grad_head_pva.setZero(dims, 3);
  result.grad_tail_pva.setZero(dims, 3);
  result.grad_durations = total_grad_t_direct;

  const Eigen::MatrixXd adj_grad = McoOptimizer::solveBandedSystemAdj(
      factorized_m, total_grad_c, kQuinticCoeffCount, kQuinticCoeffCount);

  if (adj_grad.rows() == 0) {
    return result;
  }

  for (int i = 0; i < n_segments - 1; ++i) {
    const int row_waypoint = s + i * (2 * s) + (s - 1);
    result.grad_waypoints.row(i) = adj_grad.row(row_waypoint);
  }

  result.grad_head_pva.col(0) = adj_grad.row(0).transpose();
  result.grad_head_pva.col(1) = adj_grad.row(1).transpose();
  result.grad_head_pva.col(2) = adj_grad.row(2).transpose();

  result.grad_tail_pva.col(0) = adj_grad.row(kQuinticCoeffCount * n_segments - 3).transpose();
  result.grad_tail_pva.col(1) = adj_grad.row(kQuinticCoeffCount * n_segments - 2).transpose();
  result.grad_tail_pva.col(2) = adj_grad.row(kQuinticCoeffCount * n_segments - 1).transpose();

  for (int i = 0; i < n_segments - 1; ++i) {
    const double t1 = durations(i);
    const double t2 = t1 * t1;
    const double t3 = t2 * t1;
    const double t4 = t3 * t1;

    const auto c1 = coeffs.row(i * kQuinticCoeffCount + 1);
    const auto c2 = coeffs.row(i * kQuinticCoeffCount + 2);
    const auto c3 = coeffs.row(i * kQuinticCoeffCount + 3);
    const auto c4 = coeffs.row(i * kQuinticCoeffCount + 4);
    const auto c5 = coeffs.row(i * kQuinticCoeffCount + 5);

    const Eigen::RowVectorXd neg_vel =
        -(c1 + 2.0 * t1 * c2 + 3.0 * t2 * c3 + 4.0 * t3 * c4 + 5.0 * t4 * c5);
    const Eigen::RowVectorXd neg_acc =
        -(2.0 * c2 + 6.0 * t1 * c3 + 12.0 * t2 * c4 + 20.0 * t3 * c5);
    const Eigen::RowVectorXd neg_jerk =
        -(6.0 * c3 + 24.0 * t1 * c4 + 60.0 * t2 * c5);
    const Eigen::RowVectorXd neg_snap =
        -(24.0 * c4 + 120.0 * t1 * c5);
    const Eigen::RowVectorXd neg_crackle = -120.0 * c5;

    const int row = 2 * s * i + s;
    result.grad_durations(i) +=
        neg_snap.dot(adj_grad.row(row + 0)) +
        neg_crackle.dot(adj_grad.row(row + 1)) +
        neg_vel.dot(adj_grad.row(row + 2)) +
        neg_vel.dot(adj_grad.row(row + 3)) +
        neg_acc.dot(adj_grad.row(row + 4)) +
        neg_jerk.dot(adj_grad.row(row + 5));
  }

  if (n_segments > 0) {
    const int i = n_segments - 1;
    const double t1 = durations(i);
    const double t2 = t1 * t1;
    const double t3 = t2 * t1;
    const double t4 = t3 * t1;

    const auto c1 = coeffs.row(i * kQuinticCoeffCount + 1);
    const auto c2 = coeffs.row(i * kQuinticCoeffCount + 2);
    const auto c3 = coeffs.row(i * kQuinticCoeffCount + 3);
    const auto c4 = coeffs.row(i * kQuinticCoeffCount + 4);
    const auto c5 = coeffs.row(i * kQuinticCoeffCount + 5);

    const Eigen::RowVectorXd neg_vel =
        -(c1 + 2.0 * t1 * c2 + 3.0 * t2 * c3 + 4.0 * t3 * c4 + 5.0 * t4 * c5);
    const Eigen::RowVectorXd neg_acc =
        -(2.0 * c2 + 6.0 * t1 * c3 + 12.0 * t2 * c4 + 20.0 * t3 * c5);
    const Eigen::RowVectorXd neg_jerk =
        -(6.0 * c3 + 24.0 * t1 * c4 + 60.0 * t2 * c5);

    result.grad_durations(i) +=
        neg_vel.dot(adj_grad.row(kQuinticCoeffCount * n_segments - 3)) +
        neg_acc.dot(adj_grad.row(kQuinticCoeffCount * n_segments - 2)) +
        neg_jerk.dot(adj_grad.row(kQuinticCoeffCount * n_segments - 1));
  }

  return result;
}

double accumulateSampledCostAndGradient(const AckermannBilayerParameters& params,
                                        const BilayerContext& context,
                                        Eigen::MatrixXd& grad_c_geom,
                                        Eigen::VectorXd& grad_dS_direct_geom,
                                        Eigen::MatrixXd& grad_c_scale_dummy,
                                        Eigen::VectorXd& grad_dT_direct)
{
  const int n_segments = params.n_segments;
  const int samples_per_segment = std::max(2, params.kappa);
  const double eta = static_cast<double>(params.gear_sign);
  const double L = params.wheel_base;
  const std::vector<Eigen::Vector2d> footprint =
      params.footprint_points.empty() ? std::vector<Eigen::Vector2d>{Eigen::Vector2d::Zero()}
                                      : params.footprint_points;

  grad_c_geom.setZero(context.geom_coeffs.rows(), 2);
  grad_dS_direct_geom.setZero(n_segments);
  grad_c_scale_dummy.setZero(context.scale_coeffs_dummy.rows(), 2);
  grad_dT_direct.setZero(n_segments);

  double total_cost = 0.0;
  Eigen::Matrix2d B;
  B << 0.0, -1.0,
       1.0,  0.0;

  double base_s = 0.0;
  for (int seg = 0; seg < n_segments; ++seg) {
    const double step = context.dT(seg) / static_cast<double>(samples_per_segment);

    for (int j = 0; j <= samples_per_segment; ++j) {
      const double alpha = static_cast<double>(j) / static_cast<double>(samples_per_segment);
      const double sample_weight = (j == 0 || j == samples_per_segment) ? 0.5 : 1.0;
      const double local_t = alpha * context.dT(seg);

      const Poly1DEval scale_eval = evalScaleSegment(context.scale_coeffs_dummy, seg, local_t);
      const double local_s = scale_eval.p - base_s;
      if (local_s < -1.0e-6 || local_s > context.dS(seg) + 1.0e-6) {
        return kHugeCost;
      }

      const Poly2DEval geom_eval = evalGeometrySegment(context.geom_coeffs, seg, local_s);
      const Eigen::Vector2d& pos = geom_eval.p;
      const Eigen::Vector2d& q = geom_eval.d1;
      const Eigen::Vector2d& a = geom_eval.d2;
      const Eigen::Vector2d& jg = geom_eval.d3;
      const Eigen::Vector2d& n = geom_eval.d4;

      const double u = scale_eval.d1;
      const double r = scale_eval.d2;
      const double w = scale_eval.d3;
      const double x = scale_eval.d4;

      const double m = q.norm();
      if (!std::isfinite(m) || m <= kMinSafeNorm) {
        return kHugeCost;
      }

      const Eigen::Vector2d Bq = B * q;
      const double C = a.dot(Bq);
      const double D = a.dot(q);
      const double A = jg.dot(Bq);
      const double N = A * m * m * m - 3.0 * C * D * m;
      const double M = std::pow(m, 6) + (L * C) * (L * C);
      if (!std::isfinite(M) || M <= kMinSafeNorm) {
        return kHugeCost;
      }

      const double kappa = eta * C / (m * m * m);
      const double phi = std::atan(L * kappa);
      const double phi_dot = eta * L * N / M * u;
      const double a_lon = r * m + u * u * D / m;
      const double a_lat = eta * u * u * C / m;
      const Eigen::Vector2d sigma_jerk = q * w + 3.0 * a * u * r + jg * u * u * u;

      Eigen::Vector2d grad_pos = Eigen::Vector2d::Zero();
      Eigen::Vector2d grad_q = Eigen::Vector2d::Zero();
      Eigen::Vector2d grad_a = Eigen::Vector2d::Zero();
      Eigen::Vector2d grad_j = Eigen::Vector2d::Zero();
      double grad_u = 0.0;
      double grad_r = 0.0;
      double grad_w = 0.0;
      double ell = 0.0;

      if (params.w_jerk > 0.0) {
        ell += params.w_jerk * sigma_jerk.squaredNorm();
        grad_q += params.w_jerk * (2.0 * w * sigma_jerk);
        grad_a += params.w_jerk * (6.0 * u * r * sigma_jerk);
        grad_j += params.w_jerk * (2.0 * u * u * u * sigma_jerk);
        grad_u += params.w_jerk * (2.0 * sigma_jerk.dot(3.0 * a * r + 3.0 * u * u * jg));
        grad_r += params.w_jerk * (6.0 * u * sigma_jerk.dot(a));
        grad_w += params.w_jerk * (2.0 * sigma_jerk.dot(q));
      }

      if (params.w_pseudo_speed > 0.0) {
        double penalty = 0.0;
        double dpenalty = 0.0;
        positiveCubicPenalty(params.min_pseudo_speed * params.min_pseudo_speed - q.squaredNorm(),
                             penalty, dpenalty);
        ell += params.w_pseudo_speed * penalty;
        grad_q += params.w_pseudo_speed * dpenalty * (-2.0 * q);
      }

      if (params.w_speed > 0.0) {
        double penalty = 0.0;
        double dpenalty = 0.0;
        positiveCubicPenalty(u * u * q.squaredNorm() - params.v_max * params.v_max,
                             penalty, dpenalty);
        ell += params.w_speed * penalty;
        grad_q += params.w_speed * dpenalty * (2.0 * u * u * q);
        grad_u += params.w_speed * dpenalty * (2.0 * u * q.squaredNorm());
      }

      if (params.w_speed > 0.0) {
        double penalty = 0.0;
        double dpenalty = 0.0;
        positiveCubicPenalty(-u, penalty, dpenalty);
        ell += params.w_speed * 0.1 * penalty;
        grad_u += params.w_speed * 0.1 * dpenalty * (-1.0);
      }

      const Eigen::Vector2d dC_dq = -B * a;
      const Eigen::Vector2d dC_da = B * q;
      const Eigen::Vector2d dm_dq = q / m;
      const Eigen::Vector2d dkappa_dq =
          eta * (dC_dq / std::pow(m, 3) - 3.0 * C * q / std::pow(m, 5));
      const Eigen::Vector2d dkappa_da = eta * dC_da / std::pow(m, 3);

      if (params.w_steering > 0.0) {
        double penalty = 0.0;
        double dpenalty = 0.0;
        positiveCubicPenalty(phi * phi - params.phi_max * params.phi_max, penalty, dpenalty);
        if (penalty > 0.0) {
          const double dphi_dkappa = L / (1.0 + (L * kappa) * (L * kappa));
          const double scale = params.w_steering * dpenalty * 2.0 * phi * dphi_dkappa;
          ell += params.w_steering * penalty;
          grad_q += scale * dkappa_dq;
          grad_a += scale * dkappa_da;
        }
      }

      if (params.w_steering_rate > 0.0) {
        const Eigen::Vector2d dA_dq = -B * jg;
        const Eigen::Vector2d dA_dj = B * q;
        const Eigen::Vector2d dD_dq = a;
        const Eigen::Vector2d dD_da = q;

        const Eigen::Vector2d dN_dq =
            std::pow(m, 3) * dA_dq +
            3.0 * A * m * q -
            3.0 * (m * D * dC_dq + m * C * dD_dq + C * D * dm_dq);
        const Eigen::Vector2d dN_da =
            -3.0 * (m * D * dC_da + m * C * dD_da);
        const Eigen::Vector2d dN_dj = std::pow(m, 3) * dA_dj;
        const Eigen::Vector2d dM_dq =
            6.0 * std::pow(m, 4) * q + 2.0 * L * L * C * dC_dq;
        const Eigen::Vector2d dM_da = 2.0 * L * L * C * dC_da;

        const Eigen::Vector2d dphi_dot_dq =
            eta * L * u * (M * dN_dq - N * dM_dq) / (M * M);
        const Eigen::Vector2d dphi_dot_da =
            eta * L * u * (M * dN_da - N * dM_da) / (M * M);
        const Eigen::Vector2d dphi_dot_dj = eta * L * u * dN_dj / M;
        const double dphi_dot_du = eta * L * N / M;

        double penalty = 0.0;
        double dpenalty = 0.0;
        positiveCubicPenalty(phi_dot * phi_dot - params.phi_dot_max * params.phi_dot_max,
                             penalty, dpenalty);
        ell += params.w_steering_rate * penalty;
        const double scale = params.w_steering_rate * dpenalty * 2.0 * phi_dot;
        grad_q += scale * dphi_dot_dq;
        grad_a += scale * dphi_dot_da;
        grad_j += scale * dphi_dot_dj;
        grad_u += scale * dphi_dot_du;
      }

      if (params.w_longitudinal_acc > 0.0) {
        const Eigen::Vector2d da_lon_dq =
            r * q / m +
            u * u * (a / m - D * q / std::pow(m, 3));
        const Eigen::Vector2d da_lon_da = u * u * q / m;
        const double da_lon_du = 2.0 * u * D / m;
        const double da_lon_dr = m;

        double penalty = 0.0;
        double dpenalty = 0.0;
        positiveCubicPenalty(a_lon * a_lon - params.a_lon_max * params.a_lon_max,
                             penalty, dpenalty);
        ell += params.w_longitudinal_acc * penalty;
        const double scale = params.w_longitudinal_acc * dpenalty * 2.0 * a_lon;
        grad_q += scale * da_lon_dq;
        grad_a += scale * da_lon_da;
        grad_u += scale * da_lon_du;
        grad_r += scale * da_lon_dr;
      }

      if (params.w_lateral_acc > 0.0) {
        const Eigen::Vector2d da_lat_dq =
            eta * u * u * (dC_dq / m - C * q / std::pow(m, 3));
        const Eigen::Vector2d da_lat_da = eta * u * u * dC_da / m;
        const double da_lat_du = 2.0 * eta * u * C / m;

        double penalty = 0.0;
        double dpenalty = 0.0;
        positiveCubicPenalty(a_lat * a_lat - params.a_lat_max * params.a_lat_max,
                             penalty, dpenalty);
        ell += params.w_lateral_acc * penalty;
        const double scale = params.w_lateral_acc * dpenalty * 2.0 * a_lat;
        grad_q += scale * da_lat_dq;
        grad_a += scale * da_lat_da;
        grad_u += scale * da_lat_du;
      }

      if (params.w_obstacle > 0.0 && params.esdf_query) {
        for (const Eigen::Vector2d& local_pt : footprint) {
          Eigen::Matrix2d Sr;
          Sr << local_pt.x(), -local_pt.y(),
                local_pt.y(),  local_pt.x();
          const Eigen::Vector2d rotated_pt = eta / m * (Sr * q);
          const Eigen::Vector2d world_pt = pos + rotated_pt;

          double dist = 0.0;
          Eigen::Vector2d grad_esdf = Eigen::Vector2d::Zero();
          if (params.esdf_query(world_pt, dist, grad_esdf)) {
            double penalty = 0.0;
            double dpenalty = 0.0;
            positiveCubicPenalty(params.safe_distance - dist, penalty, dpenalty);
            ell += params.w_obstacle * penalty;

            const Eigen::Matrix2d jac_q =
                eta * Sr / m - (rotated_pt * q.transpose()) / (m * m);
            const Eigen::Vector2d obstacle_grad = -params.w_obstacle * dpenalty * grad_esdf;
            grad_pos += obstacle_grad;
            grad_q += jac_q.transpose() * obstacle_grad;
          }
        }
      }

      const double weight_step = sample_weight * step;
      total_cost += weight_step * ell;

      const Eigen::VectorXd beta0_geom = McoOptimizer::getPolyBasis(local_s, 5, 0);
      const Eigen::VectorXd beta1_geom = McoOptimizer::getPolyBasis(local_s, 5, 1);
      const Eigen::VectorXd beta2_geom = McoOptimizer::getPolyBasis(local_s, 5, 2);
      const Eigen::VectorXd beta3_geom = McoOptimizer::getPolyBasis(local_s, 5, 3);

      grad_c_geom.block(seg * kQuinticCoeffCount, 0, kQuinticCoeffCount, 2) +=
          weight_step *
          (beta0_geom * grad_pos.transpose() +
           beta1_geom * grad_q.transpose() +
           beta2_geom * grad_a.transpose() +
           beta3_geom * grad_j.transpose());

      const double grad_s_geom =
          grad_pos.dot(q) +
          grad_q.dot(a) +
          grad_a.dot(jg) +
          grad_j.dot(n);

      if (seg > 0) {
        grad_dS_direct_geom.head(seg).array() -= weight_step * grad_s_geom;
      }

      const Eigen::VectorXd beta0_scale = McoOptimizer::getPolyBasis(local_t, 5, 0);
      const Eigen::VectorXd beta1_scale = McoOptimizer::getPolyBasis(local_t, 5, 1);
      const Eigen::VectorXd beta2_scale = McoOptimizer::getPolyBasis(local_t, 5, 2);
      const Eigen::VectorXd beta3_scale = McoOptimizer::getPolyBasis(local_t, 5, 3);

      grad_c_scale_dummy.block(seg * kQuinticCoeffCount, 0, kQuinticCoeffCount, 1) +=
          weight_step *
          (beta0_scale * grad_s_geom +
           beta1_scale * grad_u +
           beta2_scale * grad_r +
           beta3_scale * grad_w);

      grad_dT_direct(seg) += sample_weight *
          (ell / static_cast<double>(samples_per_segment) +
           step * alpha * (grad_s_geom * u + grad_u * r + grad_r * w + grad_w * x));
    }

    base_s += context.dS(seg);
  }

  return total_cost;
}

bool buildOutputTrajectory(const AckermannBilayerParameters& params,
                           const Eigen::VectorXd& x_vec,
                           AckermannBilayerTrajectory& trajectory)
{
  BilayerContext context;
  if (!solveBilayerTrajectory(params, x_vec, context)) {
    return false;
  }

  trajectory.geom_coeffs = context.geom_coeffs;
  trajectory.time_coeffs = context.scale_coeffs_dummy.block(0, 0, context.scale_coeffs_dummy.rows(), 1);
  trajectory.pseudo_durations = context.dS;
  trajectory.time_durations = context.dT;
  trajectory.total_duration = context.dT.sum();
  trajectory.total_pseudo_length = context.dS.sum();
  trajectory.final_cost = 0.0;
  trajectory.gear_sign = params.gear_sign;
  trajectory.waypoints = context.waypoints;
  return trajectory.isValid();
}

}  // namespace

Eigen::VectorXd AckermannBilayerOptimizer::forwardT(const Eigen::VectorXd& tau)
{
  return McoOptimizer::forwardT(tau);
}

Eigen::VectorXd AckermannBilayerOptimizer::backwardT(const Eigen::VectorXd& T)
{
  return McoOptimizer::backwardT(T);
}

Eigen::VectorXd AckermannBilayerOptimizer::backwardGradT(const Eigen::VectorXd& tau,
                                                         const Eigen::VectorXd& gradT)
{
  return McoOptimizer::backwardGradT(tau, gradT);
}

bool AckermannBilayerOptimizer::optimize(const AckermannBilayerParameters& params,
                                         AckermannBilayerTrajectory& trajectory)
{
  const int n_segments = params.n_segments;
  const int n_intermediate = std::max(0, n_segments - 1);
  const int spatial_len = params.use_vmc ? n_intermediate : 2 * n_intermediate;

  if (n_segments <= 0 ||
      params.gear_sign == 0 ||
      static_cast<int>(params.prey_points.size()) != n_segments + 1) {
    return false;
  }

  Eigen::VectorXd dT0(n_segments);
  Eigen::VectorXd dS0(n_segments);
  for (int i = 0; i < n_segments; ++i) {
    const double segment_length =
        (params.prey_points[static_cast<size_t>(i) + 1U] - params.prey_points[static_cast<size_t>(i)]).norm();
    dT0(i) = (static_cast<int>(params.initial_time_durations.size()) == n_segments)
        ? std::max(1.0e-3, params.initial_time_durations[static_cast<size_t>(i)])
        : std::max(0.2, segment_length / std::max(0.2, params.v_max * 0.7));
    dS0(i) = (static_cast<int>(params.initial_pseudo_durations.size()) == n_segments)
        ? std::max(1.0e-3, params.initial_pseudo_durations[static_cast<size_t>(i)])
        : std::max(params.boundary_tangent_norm, segment_length);
  }

  Eigen::VectorXd x(2 * n_segments + spatial_len);
  x.head(n_segments) = backwardT(dT0);
  x.segment(n_segments, n_segments) = backwardT(dS0);
  if (spatial_len > 0) {
    if (params.use_vmc) {
      x.tail(spatial_len).setOnes();
    } else {
      for (int i = 0; i < n_intermediate; ++i) {
        x(2 * n_segments + 2 * i) = params.prey_points[static_cast<size_t>(i) + 1U].x();
        x(2 * n_segments + 2 * i + 1) = params.prey_points[static_cast<size_t>(i) + 1U].y();
      }
    }
  }

  LbfgsData data{&params};
  lbfgs::lbfgs_parameter_t lbfgs_param;
  lbfgs_param.max_iterations = 180;
  lbfgs_param.g_epsilon = 1.0e-7;
  lbfgs_param.mem_size = 12;
  lbfgs_param.max_linesearch = 64;
  lbfgs_param.f_dec_coeff = 1.0e-4;
  lbfgs_param.s_curv_coeff = 0.9;

  struct CallbackData {
    LbfgsData* data{nullptr};
  } callback_data{&data};

  auto eval = [](void* instance, const Eigen::VectorXd& xin, Eigen::VectorXd& gout) -> double {
    auto* cb = reinterpret_cast<CallbackData*>(instance);
    return AckermannBilayerOptimizer::costFunctionVector(xin, &gout, cb->data);
  };

  double final_cost = 0.0;
  Eigen::VectorXd x_opt = x;
  Eigen::VectorXd g_init = Eigen::VectorXd::Zero(x_opt.size());
  (void)costFunctionVector(x_opt, &g_init, &data);
  (void)lbfgs::lbfgs_optimize(x_opt, final_cost, eval, nullptr, nullptr, &callback_data, lbfgs_param);
  final_cost = costFunctionVector(x_opt, nullptr, &data);

  if (!buildOutputTrajectory(params, x_opt, trajectory)) {
    return false;
  }
  trajectory.final_cost = final_cost;
  return std::isfinite(final_cost) && final_cost < kHugeCost;
}

double AckermannBilayerOptimizer::costFunctionVector(const Eigen::VectorXd& x_vec,
                                                     Eigen::VectorXd* g_out,
                                                     void* data)
{
  const auto* wrapper = reinterpret_cast<LbfgsData*>(data);
  const AckermannBilayerParameters& params = *wrapper->params;

  const int n_segments = params.n_segments;
  const int n_intermediate = std::max(0, n_segments - 1);
  const int spatial_len = params.use_vmc ? n_intermediate : 2 * n_intermediate;

  if (x_vec.size() != 2 * n_segments + spatial_len) {
    if (g_out) {
      g_out->setZero(x_vec.size());
    }
    return kHugeCost;
  }

  BilayerContext context;
  if (!solveBilayerTrajectory(params, x_vec, context)) {
    if (g_out) {
      g_out->setZero(x_vec.size());
    }
    return kHugeCost;
  }

  double distance_cost = 0.0;
  Eigen::MatrixXd grad_q_distance;
  buildDistanceGradient(context.waypoints, params.w_distance, grad_q_distance, distance_cost);

  Eigen::MatrixXd grad_c_geom;
  Eigen::VectorXd grad_dS_direct_geom;
  Eigen::MatrixXd grad_c_scale_dummy;
  Eigen::VectorXd grad_dT_direct;
  const double sampled_cost = accumulateSampledCostAndGradient(
      params, context,
      grad_c_geom, grad_dS_direct_geom,
      grad_c_scale_dummy, grad_dT_direct);
  if (!std::isfinite(sampled_cost) || sampled_cost >= kHugeCost) {
    if (g_out) {
      g_out->setZero(x_vec.size());
    }
    return kHugeCost;
  }

  const double time_cost = params.w_time * context.dT.sum();
  const double total_cost = sampled_cost + distance_cost + time_cost;

  if (!g_out) {
    return total_cost;
  }

  const MincoBackprop geom_backprop = backpropagateMinco(
      context.dS, grad_c_geom, grad_dS_direct_geom, context.geom_factor, context.geom_coeffs);

  const MincoBackprop scale_backprop = backpropagateMinco(
      context.dT,
      grad_c_scale_dummy,
      grad_dT_direct + params.w_time * Eigen::VectorXd::Ones(n_segments),
      context.scale_factor,
      context.scale_coeffs_dummy);

  Eigen::MatrixXd grad_q_total = geom_backprop.grad_waypoints;
  if (grad_q_distance.size() == grad_q_total.size()) {
    grad_q_total += grad_q_distance;
  }

  Eigen::VectorXd grad_dS_from_scale = Eigen::VectorXd::Zero(n_segments);
  double running = scale_backprop.grad_tail_pva(0, 0);
  grad_dS_from_scale(n_segments - 1) = running;
  for (int i = n_segments - 2; i >= 0; --i) {
    running += scale_backprop.grad_waypoints(i, 0);
    grad_dS_from_scale(i) = running;
  }

  const Eigen::VectorXd grad_dS_total = geom_backprop.grad_durations + grad_dS_from_scale;
  const Eigen::VectorXd tauT = x_vec.head(n_segments);
  const Eigen::VectorXd tauS = x_vec.segment(n_segments, n_segments);
  const Eigen::VectorXd grad_tauT = backwardGradT(tauT, scale_backprop.grad_durations);
  const Eigen::VectorXd grad_tauS = backwardGradT(tauS, grad_dS_total);
  const Eigen::VectorXd grad_spatial =
      spatial_len > 0 ? backpropagateGradQtoSpatial(grad_q_total, params) : Eigen::VectorXd();

  g_out->resize(x_vec.size());
  g_out->head(n_segments) = grad_tauT;
  g_out->segment(n_segments, n_segments) = grad_tauS;
  if (spatial_len > 0) {
    g_out->tail(spatial_len) = grad_spatial;
  }

  if (!g_out->allFinite()) {
    g_out->setZero(x_vec.size());
    return kHugeCost;
  }

  return total_cost;
}

AckermannBodyState AckermannBilayerOptimizer::evaluateBodyState(const AckermannBilayerTrajectory& trajectory,
                                                                const double wheel_base,
                                                                const double t)
{
  AckermannBodyState state;
  if (!trajectory.isValid()) {
    return state;
  }

  const SegmentSample time_sample = locateSegmentSample(trajectory.time_durations, t);

  double base_s = 0.0;
  for (int i = 0; i < time_sample.index; ++i) {
    base_s += trajectory.pseudo_durations(i);
  }

  const Poly1DEval scale_eval = evalScaleSegment(trajectory.time_coeffs, time_sample.index, time_sample.local_time);
  const double local_s = scale_eval.p - base_s;
  const Poly2DEval geom_eval = evalGeometrySegment(trajectory.geom_coeffs, time_sample.index, local_s);

  Eigen::Matrix2d B;
  B << 0.0, -1.0,
       1.0,  0.0;

  const double m = std::max(kMinSafeNorm, geom_eval.d1.norm());
  const double eta = static_cast<double>(trajectory.gear_sign);
  const double C = geom_eval.d2.dot(B * geom_eval.d1);
  const double D = geom_eval.d2.dot(geom_eval.d1);
  const double A = geom_eval.d3.dot(B * geom_eval.d1);
  const double N = A * m * m * m - 3.0 * C * D * m;
  const double M = std::pow(m, 6) + std::pow(wheel_base * C, 2);

  state.position = geom_eval.p;
  state.yaw = std::atan2(eta * geom_eval.d1.y(), eta * geom_eval.d1.x());
  state.speed = eta * scale_eval.d1 * m;
  state.curvature = eta * C / std::pow(m, 3);
  state.steering = std::atan(wheel_base * state.curvature);
  state.steering_rate = eta * wheel_base * N / M * scale_eval.d1;
  state.longitudinal_acc = scale_eval.d2 * m + scale_eval.d1 * scale_eval.d1 * D / m;
  state.lateral_acc = eta * scale_eval.d1 * scale_eval.d1 * C / m;
  return state;
}

std::vector<Eigen::Vector2d> AckermannBilayerOptimizer::samplePositions(
    const AckermannBilayerTrajectory& trajectory,
    const double dt)
{
  std::vector<Eigen::Vector2d> points;
  if (!trajectory.isValid() || dt <= 0.0) {
    return points;
  }

  const double total = trajectory.total_duration;
  for (double t = 0.0; t < total; t += dt) {
    points.push_back(evaluateBodyState(trajectory, 1.0, t).position);
  }
  points.push_back(evaluateBodyState(trajectory, 1.0, total).position);
  return points;
}

}  // namespace mocha
