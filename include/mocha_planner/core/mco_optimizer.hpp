#pragma once

/**
 * @file mco_optimizer.hpp
 * @brief MCO trajectory optimiser declarations
 * @details No ROS dependencies.  Uses L-BFGS for unconstrained optimisation.
 */

#include <vector>
#include <cmath>
#include <Eigen/Core>
#include "lbfgs.hpp"
#include "mocha_planner/core/types.hpp"

namespace mocha {

/**
 * @class McoOptimizer
 * @brief Static-method class implementing the MINCO trajectory optimiser
 */
class McoOptimizer {
public:

  /// Data block forwarded through the L-BFGS callback
  struct LbfgsData {
    const McoParameters* params;          ///< parameter pointer
    Eigen::MatrixXd M_banded;             ///< banded matrix storage
  };

  /// Main entry point: L-BFGS over x=[tau; s]
  static bool optimize(const McoParameters& params, McoTrajectory& trajectory);

  /// Cost function (vector form) and gradient
  static double costFunctionVector(const Eigen::VectorXd& x_vec, Eigen::VectorXd* g_out, void* data);

  /// Build the MCO linear system  M * c = b
  static void buildMcoMatrix(const std::vector<Eigen::Vector2d>& waypoints,
                             const Eigen::VectorXd& T,
                             const McoParameters& params,
                             Eigen::MatrixXd& M_banded,
                             Eigen::MatrixXd& b);

  /// Solve banded linear system (in-place LU)
  static Eigen::MatrixXd solveBandedSystem(Eigen::MatrixXd& M_banded,
                                           const Eigen::MatrixXd& b,
                                           int p, int q);

  /// Solve transpose system A^T x = b using factored banded matrix
  static Eigen::MatrixXd solveBandedSystemAdj(const Eigen::MatrixXd& M_banded,
                                              const Eigen::MatrixXd& b,
                                              int p, int q);

  /// Energy cost and gradient w.r.t. coefficients and segment times
  static void calculateEnergyAndGradient(const Eigen::MatrixXd& coeffs,
                                         const Eigen::VectorXd& T,
                                         const McoParameters& params,
                                         double& energy,
                                         Eigen::MatrixXd& grad_c,
                                         Eigen::VectorXd& grad_T);

  /// Feasibility / obstacle penalty and gradient
  static void calculatePenaltyAndGradient(const Eigen::MatrixXd& coeffs,
                                          const Eigen::VectorXd& T,
                                          const McoParameters& params,
                                          double& penalty,
                                          Eigen::MatrixXd& grad_c,
                                          Eigen::VectorXd& grad_T);

  /// Back-propagate gradient through the MCO linear system
  static void backpropagateMcoGradient(const Eigen::VectorXd& T,
                                       const Eigen::MatrixXd& total_grad_c,
                                       const Eigen::VectorXd& total_grad_T_direct,
                                       const McoParameters& params,
                                       const Eigen::MatrixXd& M_banded_LU,
                                       const Eigen::MatrixXd& coeffs,
                                       Eigen::MatrixXd& grad_q,
                                       Eigen::VectorXd& grad_T_total);

  /// Gradient propagation from waypoints q to motion-camouflage variables v
  static Eigen::VectorXd backpropagateGradQtoV(const Eigen::MatrixXd& grad_q,
                                               const McoParameters& params);

  // Time <-> unconstrained variable mappings
  static Eigen::VectorXd forwardT(const Eigen::VectorXd& tau);
  static Eigen::VectorXd backwardT(const Eigen::VectorXd& T);
  static Eigen::VectorXd backwardGradT(const Eigen::VectorXd& tau, const Eigen::VectorXd& gradT);

  // Polynomial basis evaluation
  static void getPolyBases(double t, int n_order,
                           Eigen::VectorXd& B0, Eigen::VectorXd& B1,
                           Eigen::VectorXd& B2, Eigen::VectorXd& B3);
  static Eigen::VectorXd getPolyBasis(double t, int n_order, int derivative_order);
};

} // namespace mocha
