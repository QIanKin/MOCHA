#include "mocha_planner/core/mco_optimizer.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <limits>

namespace mocha {

namespace {

constexpr int kQuinticCoeffCount = 6;
constexpr double kDistanceSmoothEps = 1.0e-6;

std::vector<Eigen::Vector2d> buildVmcWaypoints(const Eigen::VectorXd& v_vec,
                                               const McoParameters& params)
{
  const int n_segments = params.n_segments;
  const int v_len = std::max(0, n_segments - 1);

  std::vector<Eigen::Vector2d> waypoints;
  waypoints.reserve(n_segments + 1);
  waypoints.push_back(params.start_waypoint);

  for (int i = 0; i < v_len; ++i) {
    const Eigen::Vector2d prey_i = params.prey_points[i + 1];
    waypoints.push_back(params.ref_point + v_vec(i) * (prey_i - params.ref_point));
  }

  waypoints.push_back(params.end_waypoint);
  return waypoints;
}

std::vector<Eigen::Vector2d> buildDirectWaypoints(const Eigen::VectorXd& q_vec,
                                                   const McoParameters& params)
{
  const int n_segments = params.n_segments;
  const int n_intermediate = std::max(0, n_segments - 1);

  std::vector<Eigen::Vector2d> waypoints;
  waypoints.reserve(n_segments + 1);
  waypoints.push_back(params.start_waypoint);

  for (int i = 0; i < n_intermediate; ++i) {
    waypoints.push_back(Eigen::Vector2d(q_vec(2 * i), q_vec(2 * i + 1)));
  }

  waypoints.push_back(params.end_waypoint);
  return waypoints;
}

}  // namespace

// ----- q <- v backprop -----
Eigen::VectorXd McoOptimizer::backpropagateGradQtoV(const Eigen::MatrixXd& grad_q,
                                                    const McoParameters& params) {
  const int n_intermediate = params.n_segments > 1 ? params.n_segments - 1 : 0;
  Eigen::VectorXd grad_v = Eigen::VectorXd::Zero(n_intermediate);

  for (int i = 0; i < n_intermediate; ++i) {
    if ((int)params.prey_points.size() <= i + 1) {
      return grad_v; // zeros
    }
    Eigen::Vector2d prey_i = params.prey_points[i + 1];
    Eigen::Vector2d dq_dv = prey_i - params.ref_point;
    double gv = 0.0;
    if (grad_q.cols() >= 2) {
      gv = grad_q(i, 0) * dq_dv.x() + grad_q(i, 1) * dq_dv.y();
    }
    grad_v(i) = gv;
  }
  return grad_v;
}

// ----- optimize main -----
bool McoOptimizer::optimize(const McoParameters& params, McoTrajectory& trajectory)
{
  const int n_segments = params.n_segments;
  const int n_intermediate = std::max(0, n_segments - 1);
  const bool vmc = params.use_vmc;
  // VMC: spatial_len = n_intermediate (v params)
  // Direct: spatial_len = 2 * n_intermediate (qx,qy pairs)
  const int spatial_len = vmc ? n_intermediate : 2 * n_intermediate;

  if (n_segments <= 0 || static_cast<int>(params.prey_points.size()) < n_segments + 1) {
    return false;
  }

  // 1) initialise segment times
  Eigen::VectorXd T_initial(n_segments);
  if ((int)params.initial_segment_times.size() == n_segments) {
    for (int i = 0; i < n_segments; ++i) {
      T_initial(i) = std::max(1e-6, params.initial_segment_times[i]);
    }
  } else {
    const double reasonable_speed = std::max(1e-6, std::min(params.v_max * 0.5, 2.0));
    for (int i = 0; i < n_segments; ++i) {
      const double seg_length = (params.prey_points[i + 1] - params.prey_points[i]).norm();
      const double T_cruise = seg_length / reasonable_speed;
      const double T_phys = 2.0 * seg_length / std::max(1e-6, params.v_max);
      T_initial(i) = std::max(0.3, std::max(T_cruise, T_phys));
    }
  }
  Eigen::VectorXd tau_initial = backwardT(T_initial);

  // 2) spatial initial values
  Eigen::VectorXd spatial_initial(spatial_len);
  if (vmc) {
    // v=1 means on prey path
    spatial_initial.setOnes();
  } else {
    // directly use prey_points xy as initial guess
    for (int i = 0; i < n_intermediate; ++i) {
      spatial_initial(2 * i)     = params.prey_points[i + 1].x();
      spatial_initial(2 * i + 1) = params.prey_points[i + 1].y();
    }
  }

  // x=[tau; spatial]
  Eigen::VectorXd x(n_segments + spatial_len);
  x.head(n_segments) = tau_initial;
  if (spatial_len > 0) x.tail(spatial_len) = spatial_initial;

  // L-BFGS
  LbfgsData data_wrapper{&params, Eigen::MatrixXd()};
  lbfgs::lbfgs_parameter_t lbfgs_param;
  lbfgs_param.max_iterations = std::max(1, params.lbfgs_max_iterations);
  lbfgs_param.g_epsilon = 1e-6;
  lbfgs_param.mem_size=10;
  lbfgs_param.max_linesearch=64;
  lbfgs_param.f_dec_coeff=1e-4;
  lbfgs_param.s_curv_coeff=0.9;

  double final_cost = 0.0;
  auto fg = [&](const Eigen::VectorXd& xvec, Eigen::VectorXd& gvec) -> double {
    return McoOptimizer::costFunctionVector(xvec, &gvec, &data_wrapper);
  };
  struct CallbackData { decltype(fg)* fn; } cbd{ &fg };
  auto eval = [](void* instance, const Eigen::VectorXd& x_in, Eigen::VectorXd& g_out) -> double {
    auto* data = reinterpret_cast<CallbackData*>(instance);
    return (*(data->fn))(x_in, g_out);
  };

  Eigen::VectorXd xvec = x;
  {
    Eigen::VectorXd g_init = Eigen::VectorXd::Zero(xvec.size());
    (void)McoOptimizer::costFunctionVector(xvec, &g_init, &data_wrapper);
  }
  int ret = lbfgs::lbfgs_optimize(xvec, final_cost, eval, nullptr, nullptr, &cbd, lbfgs_param);
  (void)ret;
  x = xvec;

  // decode
  Eigen::VectorXd tau_opt = x.head(n_segments);
  trajectory.T = forwardT(tau_opt);
  trajectory.final_cost = final_cost;

  std::vector<Eigen::Vector2d> final_wps;
  if (spatial_len > 0) {
    Eigen::VectorXd spatial_opt = x.tail(spatial_len);
    final_wps = vmc ? buildVmcWaypoints(spatial_opt, params)
                    : buildDirectWaypoints(spatial_opt, params);
  } else {
    final_wps = {params.start_waypoint, params.end_waypoint};
  }

  Eigen::MatrixXd M_banded_final, b_final;
  buildMcoMatrix(final_wps, trajectory.T, params, M_banded_final, b_final);
  trajectory.coeffs = solveBandedSystem(M_banded_final, b_final, params.n_coeffs, params.n_coeffs);
  trajectory.total_duration = trajectory.T.sum();
  return trajectory.coeffs.allFinite() && trajectory.T.allFinite() && std::isfinite(final_cost);
}

// ----- objective -----
double McoOptimizer::costFunctionVector(const Eigen::VectorXd& x_vec, Eigen::VectorXd* g_out, void* data)
{
  LbfgsData* d = reinterpret_cast<LbfgsData*>(data);
  const McoParameters& params = *d->params;
  const int n_segments = params.n_segments;
  const int n_intermediate = std::max(0, n_segments - 1);
  const bool vmc = params.use_vmc;
  const int spatial_len = vmc ? n_intermediate : 2 * n_intermediate;
  const int tau_len = n_segments;

  Eigen::VectorXd tau = x_vec.head(tau_len);
  Eigen::VectorXd spatial = spatial_len > 0 ? x_vec.tail(spatial_len) : Eigen::VectorXd();

  Eigen::VectorXd T = forwardT(tau);

  if (!tau.allFinite() || !T.allFinite() || (spatial_len > 0 && !spatial.allFinite()) ||
      static_cast<int>(params.prey_points.size()) < n_segments + 1) {
    if (g_out) g_out->setZero(x_vec.size());
    return 1e10;
  }

  // build waypoints
  std::vector<Eigen::Vector2d> waypoints;
  if (spatial_len > 0) {
    waypoints = vmc ? buildVmcWaypoints(spatial, params)
                    : buildDirectWaypoints(spatial, params);
  } else {
    waypoints = {params.start_waypoint, params.end_waypoint};
  }
  for (const auto& q_i : waypoints) {
    if (!q_i.allFinite()) {
      if (g_out) g_out->setZero(x_vec.size());
      return 1e10;
    }
  }

  // MCO linear system
  Eigen::MatrixXd b;
  buildMcoMatrix(waypoints, T, params, d->M_banded, b);
  Eigen::MatrixXd M_banded_factor = d->M_banded;
  Eigen::MatrixXd coeffs = solveBandedSystem(M_banded_factor, b, params.n_coeffs, params.n_coeffs);
  if (coeffs.hasNaN()) {
    if (g_out) g_out->setZero(x_vec.size());
    return 1e9;
  }

  // cost terms
  double energy_cost = 0.0, penalty_cost = 0.0;
  Eigen::MatrixXd grad_E_c, grad_P_c;
  Eigen::VectorXd grad_E_T, grad_P_T;
  calculateEnergyAndGradient(coeffs, T, params, energy_cost, grad_E_c, grad_E_T);
  calculatePenaltyAndGradient(coeffs, T, params, penalty_cost, grad_P_c, grad_P_T);
  const double time_cost = T.sum();

  // ---- waypoint distance cost ----
  double distance_cost = 0.0;
  Eigen::MatrixXd grad_q_distance;
  const int M = static_cast<int>(waypoints.size()) - 1;
  if (n_intermediate > 0 && params.w_distance > 0.0 && M >= 1) {
    grad_q_distance.setZero(n_intermediate, params.dims);
    for (int i = 0; i < M; ++i) {
      const Eigen::Vector2d delta = waypoints[i + 1] - waypoints[i];
      const double smoothed_distance = std::sqrt(delta.squaredNorm() + kDistanceSmoothEps);
      const Eigen::Vector2d grad = delta / smoothed_distance;

      distance_cost += smoothed_distance;
      if (i > 0) {
        grad_q_distance.row(i - 1) -= grad.transpose();
      }
      if (i < n_intermediate) {
        grad_q_distance.row(i) += grad.transpose();
      }
    }
    distance_cost *= params.w_distance;
    grad_q_distance *= params.w_distance;
  }

  const double total_cost = params.w_energy * energy_cost + penalty_cost
                          + params.w_time * time_cost + distance_cost;

  // analytic gradient
  if (g_out) {
    g_out->resize(tau_len + spatial_len);

    Eigen::MatrixXd total_grad_c = params.w_energy * grad_E_c + grad_P_c;
    Eigen::VectorXd total_grad_T_direct =
        params.w_energy * grad_E_T + grad_P_T +
        params.w_time * Eigen::VectorXd::Ones(n_segments);

    Eigen::MatrixXd grad_q;
    Eigen::VectorXd grad_T_total;
    backpropagateMcoGradient(T, total_grad_c, total_grad_T_direct, params, M_banded_factor, coeffs,
                             grad_q, grad_T_total);

    if (grad_q_distance.size() > 0 && grad_q.rows() == grad_q_distance.rows()) {
      grad_q += grad_q_distance;
    }

    Eigen::VectorXd grad_tau = backwardGradT(tau, grad_T_total);
    g_out->head(tau_len) = grad_tau;

    if (spatial_len > 0) {
      if (vmc) {
        // VMC: grad_q -> grad_v via chain rule
        Eigen::VectorXd grad_v = backpropagateGradQtoV(grad_q, params);
        g_out->tail(spatial_len) = grad_v;
      } else {
        // Direct MINCO: grad_q maps straight to [gqx0,gqy0,gqx1,gqy1,...]
        Eigen::VectorXd grad_spatial(spatial_len);
        for (int i = 0; i < n_intermediate; ++i) {
          grad_spatial(2 * i)     = grad_q(i, 0);
          grad_spatial(2 * i + 1) = grad_q(i, 1);
        }
        g_out->tail(spatial_len) = grad_spatial;
      }
    }
  }
  return total_cost;
}

// ----- time reparameterization -----
Eigen::VectorXd McoOptimizer::forwardT(const Eigen::VectorXd& tau) {
  Eigen::VectorXd T(tau.size());
  for (int i = 0; i < tau.size(); ++i) {
    if (tau(i) > 0) {
      T(i) = (0.5 * tau(i) + 1.0) * tau(i) + 1.0;
    } else {
      T(i) = 1.0 / ((0.5 * tau(i) - 1.0) * tau(i) + 1.0);
    }
  }
  return T;
}
Eigen::VectorXd McoOptimizer::backwardT(const Eigen::VectorXd& T) {
  Eigen::VectorXd tau(T.size());
  for (int i = 0; i < T.size(); ++i) {
    double t_val = T(i);
    if (t_val > 1.0) {
      tau(i) = std::sqrt(2.0 * t_val - 1.0) - 1.0;
    } else {
      if (t_val < 1e-12) t_val = 1e-12;
      tau(i) = 1.0 - std::sqrt(2.0 / t_val - 1.0);
    }
  }
  return tau;
}
Eigen::VectorXd McoOptimizer::backwardGradT(const Eigen::VectorXd& tau, const Eigen::VectorXd& gradT) {
  Eigen::VectorXd gradTau(tau.size());
  for (int i = 0; i < tau.size(); ++i) {
    double dT_dtau = 0.0;
    if (tau(i) > 0) dT_dtau = tau(i) + 1.0;
    else {
      const double den = (0.5 * tau(i) - 1.0) * tau(i) + 1.0;
      dT_dtau = -(tau(i) - 1.0) / (den * den);
    }
    gradTau(i) = gradT(i) * dT_dtau;
  }
  return gradTau;
}

// ----- poly bases -----
void McoOptimizer::getPolyBases(double t, int n_order,
                                Eigen::VectorXd& B0, Eigen::VectorXd& B1,
                                Eigen::VectorXd& B2, Eigen::VectorXd& B3) {
  int n_coeffs = n_order + 1;
  B0.setZero(n_coeffs);
  B1.setZero(n_coeffs);
  B2.setZero(n_coeffs);
  B3.setZero(n_coeffs);
  for (int i = 0; i <= n_order; ++i) B0(i) = std::pow(t, i);
  for (int i = 1; i <= n_order; ++i) B1(i) = i * std::pow(t, i - 1);
  for (int i = 2; i <= n_order; ++i) B2(i) = i * (i - 1) * std::pow(t, i - 2);
  for (int i = 3; i <= n_order; ++i) B3(i) = i * (i - 1) * (i - 2) * std::pow(t, i - 3);
}
Eigen::VectorXd McoOptimizer::getPolyBasis(double t, int n_order, int derivative_order) {
  int n_coeffs = n_order + 1;
  Eigen::VectorXd basis = Eigen::VectorXd::Zero(n_coeffs);
  if (derivative_order > n_order) return basis;
  for (int i = derivative_order; i <= n_order; ++i) {
    double coeff = 1.0;
    for (int j = 0; j < derivative_order; ++j) coeff *= (i - j);
    basis(i) = coeff * std::pow(t, i - derivative_order);
  }
  return basis;
}

// ----- build/solve MCO linear system -----
void McoOptimizer::buildMcoMatrix(const std::vector<Eigen::Vector2d>& waypoints,
                                  const Eigen::VectorXd& T,
                                  const McoParameters& params,
                                  Eigen::MatrixXd& M_banded,
                                  Eigen::MatrixXd& b) {
  const int n_segments = params.n_segments;
  const int n_coeffs = params.n_coeffs;
  const int total_vars = n_segments * n_coeffs;
  const int s = params.s;
  const int dims = params.dims;

  if (n_coeffs != kQuinticCoeffCount || static_cast<int>(waypoints.size()) != n_segments + 1 ||
      T.size() != n_segments) {
    M_banded.resize(0, 0);
    b.resize(0, dims);
    return;
  }

  const int p = n_coeffs, q = n_coeffs;
  M_banded.setZero(p + q + 1, total_vars);
  b.setZero(total_vars, dims);

  auto set_val = [&](int r, int c, double val) {
    if (val != 0.0) {
      M_banded(r - c + q, c) = val;
    }
  };

  int r_base = 0;
  // head constraints
  for (int k = 0; k < s; ++k) {
    Eigen::VectorXd basis = getPolyBasis(0.0, params.n_order, k);
    for (int j = 0; j < n_coeffs; ++j) set_val(r_base + k, j, basis(j));
  }
  b.row(0) = params.start_waypoint;
  b.row(1) = params.start_vel;
  b.row(2) = params.start_acc;
  r_base += s;

  // interior constraints
  for (int i = 0; i < n_segments - 1; ++i) {
    int c_base = i * n_coeffs;
    double t_i = T(i);

    for (int k = s; k < 2 * s - 1; ++k) {
      Eigen::VectorXd basis1 = getPolyBasis(t_i, params.n_order, k);
      Eigen::VectorXd basis2 = getPolyBasis(0.0, params.n_order, k);
      for (int j = 0; j < n_coeffs; ++j) {
        set_val(r_base, c_base + j, basis1(j));
        set_val(r_base, c_base + n_coeffs + j, -basis2(j));
      }
      r_base++;
    }
    // waypoint position
    Eigen::VectorXd basis_p = getPolyBasis(t_i, params.n_order, 0);
    for (int j = 0; j < n_coeffs; ++j) set_val(r_base, c_base + j, basis_p(j));
    b.row(r_base) = waypoints[i + 1];
    r_base++;

    // lower-order continuity
    for (int k = 0; k < s; ++k) {
      Eigen::VectorXd basis1 = getPolyBasis(t_i, params.n_order, k);
      Eigen::VectorXd basis2 = getPolyBasis(0.0, params.n_order, k);
      for (int j = 0; j < n_coeffs; ++j) {
        set_val(r_base, c_base + j, basis1(j));
        set_val(r_base, c_base + n_coeffs + j, -basis2(j));
      }
      r_base++;
    }
  }

  // tail constraints
  int c_base = (n_segments - 1) * n_coeffs;
  double t_M = T.tail(1)(0);
  for (int k = 0; k < s; ++k) {
    Eigen::VectorXd basis = getPolyBasis(t_M, params.n_order, k);
    for (int j = 0; j < n_coeffs; ++j) set_val(r_base + k, c_base + j, basis(j));
  }
  b.row(r_base + 0) = params.end_waypoint;
  b.row(r_base + 1) = params.end_vel;
  b.row(r_base + 2) = params.end_acc;

}

Eigen::MatrixXd McoOptimizer::solveBandedSystem(Eigen::MatrixXd& M_banded,
                                                const Eigen::MatrixXd& b,
                                                int p, int q) {
  int N = b.rows();
  if (N == 0) return Eigen::MatrixXd(0, b.cols());

  for (int k = 0; k < N - 1; ++k) {
    double pivot = M_banded(q, k);
    if (std::abs(pivot) < 1e-14) {
      return Eigen::MatrixXd::Constant(b.rows(), b.cols(), std::numeric_limits<double>::quiet_NaN());
    }
    for (int i = k + 1; i < std::min(k + p + 1, N); ++i) {
      M_banded(i - k + q, k) /= pivot;
    }
    for (int j = k + 1; j < std::min(k + q + 1, N); ++j) {
      double val = M_banded(k - j + q, j);
      if (val != 0.0) {
        for (int i = k + 1; i < std::min(k + p + 1, N); ++i) {
          M_banded(i - j + q, j) -= M_banded(i - k + q, k) * val;
        }
      }
    }
  }

  Eigen::MatrixXd y = b;
  for (int k = 0; k < N; ++k) {
    for (int i = k + 1; i < std::min(k + p + 1, N); ++i) {
      y.row(i) -= M_banded(i - k + q, k) * y.row(k);
    }
  }
  Eigen::MatrixXd x = y;
  for (int k = N - 1; k >= 0; --k) {
    x.row(k) /= M_banded(q, k);
    for (int i = std::max(0, k - q); i < k; ++i) {
      x.row(i) -= M_banded(i - k + q, k) * x.row(k);
    }
  }
  return x;
}

static Eigen::MatrixXd solveBandedAdjInternal(const Eigen::MatrixXd& M_banded,
                                              const Eigen::MatrixXd& b,
                                              int p, int q) {
  const int N = b.rows();
  if (N == 0) return Eigen::MatrixXd(0, b.cols());

  Eigen::MatrixXd x = b;

  int iM;
  for (int j = 0; j <= N - 1; ++j) {
    x.row(j) /= M_banded(q, j);
    iM = std::min(j + q, N - 1);
    for (int i = j + 1; i <= iM; ++i) {
      double aji = M_banded(j - i + q, i);
      if (aji != 0.0) {
        x.row(i) -= aji * x.row(j);
      }
    }
  }
  for (int j = N - 1; j >= 0; --j) {
    iM = std::max(0, j - p);
    for (int i = iM; i <= j - 1; ++i) {
      double aji = M_banded(j - i + q, i);
      if (aji != 0.0) {
        x.row(i) -= aji * x.row(j);
      }
    }
  }
  return x;
}

Eigen::MatrixXd McoOptimizer::solveBandedSystemAdj(const Eigen::MatrixXd& M_banded,
                                                   const Eigen::MatrixXd& b,
                                                   int p, int q) {
  return solveBandedAdjInternal(M_banded, b, p, q);
}

// ----- energy & penalty -----
void McoOptimizer::calculateEnergyAndGradient(const Eigen::MatrixXd& coeffs,
                                              const Eigen::VectorXd& T,
                                              const McoParameters& params,
                                              double& energy,
                                              Eigen::MatrixXd& grad_c,
                                              Eigen::VectorXd& grad_T) {
  energy = 0.0;
  grad_c.setZero(coeffs.rows(), coeffs.cols());
  grad_T.setZero(T.size());

  const int n_coeffs = params.n_coeffs;
  const int s = params.s;

  for (int i = 0; i < params.n_segments; ++i) {
    const auto& c3 = coeffs.row(i * n_coeffs + s);
    const auto& c4 = coeffs.row(i * n_coeffs + s + 1);
    const auto& c5 = coeffs.row(i * n_coeffs + s + 2);

    const double t1 = T(i);
    const double t2 = t1 * t1;
    const double t3 = t2 * t1;
    const double t4 = t3 * t1;
    const double t5 = t4 * t1;

    energy += 36.0 * c3.squaredNorm() * t1 +
              144.0 * c3.dot(c4) * t2 +
              192.0 * c4.squaredNorm() * t3 +
              240.0 * c3.dot(c5) * t3 +
              720.0 * c4.dot(c5) * t4 +
              720.0 * c5.squaredNorm() * t5;

    grad_c.row(i * n_coeffs + s)     = 72.0 * c3 * t1 + 144.0 * c4 * t2 + 240.0 * c5 * t3;
    grad_c.row(i * n_coeffs + s + 1) = 144.0 * c3 * t2 + 384.0 * c4 * t3 + 720.0 * c5 * t4;
    grad_c.row(i * n_coeffs + s + 2) = 240.0 * c3 * t3 + 720.0 * c4 * t4 + 1440.0 * c5 * t5;

    grad_T(i) = 36.0 * c3.squaredNorm() +
                288.0 * c3.dot(c4) * t1 +
                576.0 * c4.squaredNorm() * t2 +
                720.0 * c3.dot(c5) * t2 +
                2880.0 * c4.dot(c5) * t3 +
                3600.0 * c5.squaredNorm() * t4;
  }
}

void McoOptimizer::calculatePenaltyAndGradient(const Eigen::MatrixXd& coeffs,
                                               const Eigen::VectorXd& T,
                                               const McoParameters& params,
                                               double& penalty,
                                               Eigen::MatrixXd& grad_c,
                                               Eigen::VectorXd& grad_T) {
  penalty = 0.0;
  grad_c.setZero(coeffs.rows(), coeffs.cols());
  grad_T.setZero(T.size());

  // Prefix sums of segment durations for absolute time (dynamic obstacles)
  Eigen::VectorXd T_prefix = Eigen::VectorXd::Zero(params.n_segments);
  for (int i = 1; i < params.n_segments; ++i)
      T_prefix(i) = T_prefix(i - 1) + T(i - 1);
  // Cross-segment coupling accumulator for dynamic obstacles
  Eigen::VectorXd dyn_cross = Eigen::VectorXd::Zero(params.n_segments);

  for (int i = 0; i < params.n_segments; ++i) {
    Eigen::MatrixXd c_i = coeffs.block(i * params.n_coeffs, 0, params.n_coeffs, params.dims);
    int kappa_i = params.kappa;
    double h = T(i) / (double)kappa_i;

    for (int j = 0; j <= kappa_i; ++j) {
      double t = (double)j / (double)kappa_i * T(i);
      Eigen::VectorXd B0, B1, B2, B3;
      getPolyBases(t, params.n_order, B0, B1, B2, B3);

      Eigen::RowVector2d pos = (B0.transpose() * c_i);
      Eigen::RowVector2d vel = (B1.transpose() * c_i);
      Eigen::RowVector2d acc = (B2.transpose() * c_i);
      Eigen::RowVector2d jer = (B3.transpose() * c_i);

      double omg = (j == 0 || j == kappa_i) ? 0.5 : 1.0;

      // ---- static obstacle avoidance (ESDF): dual-threshold soft+hard ----
      if (params.esdf_query) {
        double d_margin = params.drone_radius + params.esdf_safety_margin;
        double d_soft   = d_margin + params.esdf_soft_margin;
        double dem      = params.obstacle_demarcation;
        Eigen::Vector2d pos_vec(pos(0), pos(1));
        double d_esdf = 0.0;
        Eigen::Vector2d grad_esdf = Eigen::Vector2d::Zero();
        if (params.esdf_query(pos_vec, d_esdf, grad_esdf)) {
          double hard_vio = d_margin - d_esdf;

          // --- HARD zone: d_esdf < d_margin ---
          if (hard_vio > 0.0) {
            double cost_h, dCdV;
            if (hard_vio < dem) {
              // cubic penalty for small violations
              cost_h = params.w_obstacle * hard_vio * hard_vio * hard_vio;
              dCdV   = 3.0 * params.w_obstacle * hard_vio * hard_vio;
            } else {
              // quadratic penalty for large violations (C1-continuous at dem)
              // f(dem) = w*dem^3, f'(dem) = 3*w*dem^2 (matches cubic)
              cost_h = params.w_obstacle * (3.0 * dem * hard_vio * hard_vio
                       - 3.0 * dem * dem * hard_vio + dem * dem * dem);
              dCdV   = params.w_obstacle * (6.0 * dem * hard_vio - 3.0 * dem * dem);
            }
            Eigen::RowVector2d grad_p = -dCdV * grad_esdf.transpose();

            penalty += omg * h * cost_h;
            grad_c.block(i * params.n_coeffs, 0, params.n_coeffs, params.dims) += omg * h * (B0 * grad_p);

            double alpha = (double)j / (double)kappa_i;
            double grad_T_esdf_direct = alpha * grad_p.dot(vel);
            grad_T(i) += omg * (cost_h / (double)kappa_i + h * grad_T_esdf_direct);
          }

          // --- SOFT zone: d_margin <= d_esdf < d_soft ---
          if (hard_vio <= 0.0 && d_esdf < d_soft) {
            double soft_vio = d_soft - d_esdf;
            double cost_s = params.w_obstacle_soft * soft_vio * soft_vio * soft_vio;
            Eigen::RowVector2d grad_p = -3.0 * params.w_obstacle_soft * soft_vio * soft_vio * grad_esdf.transpose();

            penalty += omg * h * cost_s;
            grad_c.block(i * params.n_coeffs, 0, params.n_coeffs, params.dims) += omg * h * (B0 * grad_p);

            double alpha = (double)j / (double)kappa_i;
            double grad_T_esdf_direct = alpha * grad_p.dot(vel);
            grad_T(i) += omg * (cost_s / (double)kappa_i + h * grad_T_esdf_direct);
          }
        }
      }

      // ---- dynamic obstacle avoidance ----
      if (!params.dynamic_obstacles.empty()) {
          const double alpha_dyn = static_cast<double>(j) / static_cast<double>(kappa_i);
          const double t_abs = params.t_plan_start + T_prefix(i) + alpha_dyn * T(i);
          const double d_base = params.drone_radius + params.dynamic_safety_margin;

          for (const auto& dob : params.dynamic_obstacles) {
              const Eigen::Vector2d p_obs =
                  dob.position + dob.velocity * (t_abs - dob.t_ref);
              const Eigen::RowVector2d delta = pos - p_obs.transpose();
              const double d_sq = delta.squaredNorm();
              const double d_safe = d_base + dob.radius;
              const double d_safe_sq = d_safe * d_safe;
              const double vio = d_safe_sq - d_sq;

              if (vio > 0.0) {
                  const double vio_sq = vio * vio;
                  const double cost_d = params.w_dynamic_obs * vio * vio_sq;
                  // grad_p = ∂cost/∂pos = -6 w vio² Δ
                  const Eigen::RowVector2d grad_p =
                      -6.0 * params.w_dynamic_obs * vio_sq * delta;
                  // relative velocity: robot minus obstacle
                  const Eigen::RowVector2d vel_rel =
                      vel - dob.velocity.transpose();

                  // accumulate cost
                  penalty += omg * h * cost_d;
                  // gradient w.r.t. polynomial coefficients
                  grad_c.block(i * params.n_coeffs, 0, params.n_coeffs, params.dims)
                      += omg * h * (B0 * grad_p);

                  // own-segment T gradient: uses vel_rel = vel_robot - v_obs
                  grad_T(i) += omg * (cost_d / static_cast<double>(kappa_i)
                               + h * alpha_dyn * grad_p.dot(vel_rel));

                  // cross-segment coupling: only -v_obs (robot velocity NOT included,
                  // because changing T_k for k<i does not affect robot's local time)
                  dyn_cross(i) += omg * h * grad_p.dot(-dob.velocity.transpose());
              }
          }
      }

      double v_pen = vel.squaredNorm() - params.v_max * params.v_max;
      if (v_pen > 0) {
        double cost_v = params.w_feasibility * std::pow(v_pen, 3);
        Eigen::RowVector2d grad_v = params.w_feasibility * 6 * std::pow(v_pen, 2) * vel;

        penalty += omg * h * cost_v;
        grad_c.block(i * params.n_coeffs, 0, params.n_coeffs, params.dims) += omg * h * (B1 * grad_v);

        double alpha = (double)j / (double)kappa_i;
        double grad_T_v_direct = alpha * grad_v.dot(acc);
        grad_T(i) += omg * (cost_v / (double)kappa_i + h * grad_T_v_direct);
      }

      double a_pen = acc.squaredNorm() - params.a_max * params.a_max;
      if (a_pen > 0) {
        double cost_a = params.w_feasibility * std::pow(a_pen, 3);
        Eigen::RowVector2d grad_a = params.w_feasibility * 6 * std::pow(a_pen, 2) * acc;

        penalty += omg * h * cost_a;
        grad_c.block(i * params.n_coeffs, 0, params.n_coeffs, params.dims) += omg * h * (B2 * grad_a);

        double alpha = (double)j / (double)kappa_i;
        double grad_T_a_direct = alpha * grad_a.dot(jer);
        grad_T(i) += omg * (cost_a / (double)kappa_i + h * grad_T_a_direct);
      }
    }
  }

  // Dynamic obstacle cross-segment time coupling — suffix sum O(n)
  // Changing T_k (k < i) shifts absolute time for all subsequent segments,
  // moving obstacle positions by v_obs * δ. Accumulated via reverse sweep.
  if (!params.dynamic_obstacles.empty()) {
      double suffix = 0.0;
      for (int i = params.n_segments - 1; i >= 0; --i) {
          grad_T(i) += suffix;
          suffix += dyn_cross(i);
      }
  }
}

// ----- backprop through MCO linear system -----
void McoOptimizer::backpropagateMcoGradient(const Eigen::VectorXd& T,
                                            const Eigen::MatrixXd& total_grad_c,
                                            const Eigen::VectorXd& total_grad_T_direct,
                                            const McoParameters& params,
                                            const Eigen::MatrixXd& M_banded_LU,
                                            const Eigen::MatrixXd& coeffs,
                                            Eigen::MatrixXd& grad_q,
                                            Eigen::VectorXd& grad_T_total) {
  int n_segments = params.n_segments;
  int s = params.s;
  int n_coeffs = params.n_coeffs;
  int dims = params.dims;

  Eigen::MatrixXd adj_grad = solveBandedSystemAdj(M_banded_LU, total_grad_c, n_coeffs, n_coeffs);

  grad_q.resize(n_segments > 1 ? n_segments - 1 : 0, dims);
  grad_T_total = total_grad_T_direct;
  for (int i = 0; i < n_segments - 1; ++i) {
    int r_waypoint = s + i * (2 * s) + (s - 1);
    grad_q.row(i) = adj_grad.row(r_waypoint);
  }

  for (int i = 0; i < n_segments - 1; ++i) {
    const double t1 = T(i);
    const double t2 = t1 * t1;
    const double t3 = t2 * t1;
    const double t4 = t3 * t1;
    const auto c1 = coeffs.row(i * n_coeffs + 1);
    const auto c2 = coeffs.row(i * n_coeffs + 2);
    const auto c3 = coeffs.row(i * n_coeffs + 3);
    const auto c4 = coeffs.row(i * n_coeffs + 4);
    const auto c5 = coeffs.row(i * n_coeffs + 5);

    const Eigen::RowVectorXd neg_vel =
        -(c1 + 2.0 * t1 * c2 + 3.0 * t2 * c3 + 4.0 * t3 * c4 + 5.0 * t4 * c5);
    const Eigen::RowVectorXd neg_acc =
        -(2.0 * c2 + 6.0 * t1 * c3 + 12.0 * t2 * c4 + 20.0 * t3 * c5);
    const Eigen::RowVectorXd neg_jerk =
        -(6.0 * c3 + 24.0 * t1 * c4 + 60.0 * t2 * c5);
    const Eigen::RowVectorXd neg_snap =
        -(24.0 * c4 + 120.0 * t1 * c5);
    const Eigen::RowVectorXd neg_crackle =
        -120.0 * c5;

    const int r = 2 * s * i + s;
    grad_T_total(i) = total_grad_T_direct(i) +
        neg_snap.dot(adj_grad.row(r + 0)) +
        neg_crackle.dot(adj_grad.row(r + 1)) +
        neg_vel.dot(adj_grad.row(r + 2)) +
        neg_vel.dot(adj_grad.row(r + 3)) +
        neg_acc.dot(adj_grad.row(r + 4)) +
        neg_jerk.dot(adj_grad.row(r + 5));
  }

  if (n_segments > 0) {
    const int i = n_segments - 1;
    const double t1 = T(i);
    const double t2 = t1 * t1;
    const double t3 = t2 * t1;
    const double t4 = t3 * t1;
    const auto c1 = coeffs.row(i * n_coeffs + 1);
    const auto c2 = coeffs.row(i * n_coeffs + 2);
    const auto c3 = coeffs.row(i * n_coeffs + 3);
    const auto c4 = coeffs.row(i * n_coeffs + 4);
    const auto c5 = coeffs.row(i * n_coeffs + 5);

    const Eigen::RowVectorXd neg_vel =
        -(c1 + 2.0 * t1 * c2 + 3.0 * t2 * c3 + 4.0 * t3 * c4 + 5.0 * t4 * c5);
    const Eigen::RowVectorXd neg_acc =
        -(2.0 * c2 + 6.0 * t1 * c3 + 12.0 * t2 * c4 + 20.0 * t3 * c5);
    const Eigen::RowVectorXd neg_jerk =
        -(6.0 * c3 + 24.0 * t1 * c4 + 60.0 * t2 * c5);

    grad_T_total(i) = total_grad_T_direct(i) +
        neg_vel.dot(adj_grad.row(n_coeffs * n_segments - 3)) +
        neg_acc.dot(adj_grad.row(n_coeffs * n_segments - 2)) +
        neg_jerk.dot(adj_grad.row(n_coeffs * n_segments - 1));
  }
}

} // namespace mocha
