#pragma once

/**
 * @file types.hpp
 * @brief MOCHA trajectory optimization core data structures
 * @details Contains parameter structs, optimization data wrappers, and trajectory results.
 *          No ROS dependencies -- only Eigen and standard C++.
 */

#include <vector>
#include <string>
#include <functional>
#include <Eigen/Core>

namespace mocha {

// ==================== ESDF query callback ====================

/**
 * @brief ESDF query function type.
 * @param pos   World position (x, y)
 * @param dist  Output: signed distance (positive = free, negative = inside obstacle)
 * @param grad  Output: normalized gradient (points away from nearest obstacle)
 * @return true if query succeeded
 */
using EsdfQueryFunc = std::function<bool(const Eigen::Vector2d& pos, double& dist, Eigen::Vector2d& grad)>;

// ==================== Dynamic obstacle model ====================

/**
 * @struct DynamicObstacle
 * @brief Linear (constant-velocity) dynamic obstacle for trajectory optimization.
 *
 * Position model:  p_obs(t) = position + velocity * (t - t_ref)
 * where t is absolute wall-clock time.
 */
struct DynamicObstacle
{
    Eigen::Vector2d position;           ///< obstacle centre at reference time (m)
    Eigen::Vector2d velocity;           ///< constant velocity (m/s)
    double t_ref = 0.0;                 ///< reference time for position (s)
    double radius = 0.3;                ///< obstacle radius (m)
};

// ==================== PRM parameters ====================

/**
 * @struct PrmParameters
 * @brief Parameters for PRM multi-homotopy path generation
 */
struct PrmParameters {
    int n_samples{200};          ///< number of random samples in local window
    int m_check{10};             ///< collision check samples per edge
    double d_safe{0.0};          ///< minimum ESDF distance for safe cells (auto = robot_radius + 0.1)
    int k_max_paths{3};          ///< maximum number of homotopy-distinct paths
    double local_range{5.0};     ///< maximum lateral sampling half-width around the guide corridor (m)
    std::vector<Eigen::Vector2d> guide_points;  ///< deterministic seed polyline, e.g. the current A* local segment
    int max_obstacle_refs{5};    ///< max number of obstacle reference points used for topology comparison
    double boundary_sample_ratio{0.3};     ///< fraction of samples biased toward obstacle boundaries
    double lateral_sample_factor{0.4};     ///< corridor half-width = clamp(factor * local_range, 2*d_safe, local_range)
    double min_guard_spacing{0.0};         ///< minimum guard separation; <=0 means derive from d_safe
    int max_raw_paths{32};                 ///< maximum number of simple graph paths explored before deduplication
    int max_path_depth{18};                ///< maximum graph depth during simple-path enumeration
    double path_length_ratio_limit{1.6};   ///< discard candidates much longer than the shortest retained path
    double smooth_spacing{0.25};           ///< temporary spacing used before topology-preserving path smoothing
    int smooth_passes{6};                  ///< number of elastic smoothing passes per candidate path
    double smooth_strength{0.5};           ///< Laplacian smoothing strength in (0, 1]
};

// ==================== MCO parameters ====================

/**
 * @struct McoParameters
 * @brief All inputs required by the MCO trajectory optimiser
 */
struct McoParameters
{
    // ---- dynamics / geometry ----
    double v_max = 2.5;               ///< maximum speed (m/s)
    double a_max = 5.0;               ///< maximum acceleration (m/s^2)
    double drone_radius = 0.20;       ///< robot radius (m)

    // ---- polynomial trajectory constants ----
    const int dims = 2;
    const int s = 3;
    const int n_order = 2 * s - 1;       // = 5
    const int n_coeffs = n_order + 1;     // = 6

    // ---- optimisation weights ----
    double w_energy = 1.0;
    double w_time = 10.0;
    double w_feasibility = 5000.0;
    double w_obstacle = 5000.0;
    double w_distance = 100.0;           ///< waypoint distance regularisation to suppress unnecessary bending

    // ---- ESDF-based collision avoidance ----
    EsdfQueryFunc esdf_query;          ///< ESDF query function (required)
    double esdf_safety_margin{0.1};    ///< additional safety margin beyond drone_radius (m)
    double w_obstacle_soft{2000.0};    ///< soft-zone penalty weight (provides far-range gradient)
    double esdf_soft_margin{0.5};      ///< soft zone extends this far beyond d_margin (m)
    double obstacle_demarcation{0.1};  ///< cubic-to-quadratic switchover depth in hard zone (m)

    // ---- dynamic obstacle avoidance ----
    std::vector<DynamicObstacle> dynamic_obstacles;   ///< moving obstacles (may be empty)
    double w_dynamic_obs = 5000.0;                    ///< dynamic obstacle penalty weight
    double dynamic_safety_margin = 0.15;              ///< additional margin beyond drone_radius + obs_radius (m)
    double t_plan_start = 0.0;                        ///< absolute time at which this plan starts (s)

    // ---- penalty sampling ----
    int kappa = 10;

    // ---- boundary conditions ----
    Eigen::Vector2d start_waypoint;
    Eigen::Vector2d end_waypoint;
    Eigen::Vector2d start_vel = Eigen::Vector2d::Zero();
    Eigen::Vector2d start_acc = Eigen::Vector2d::Zero();
    Eigen::Vector2d end_vel   = Eigen::Vector2d::Zero();
    Eigen::Vector2d end_acc   = Eigen::Vector2d::Zero();

    // ---- segment layout ----
    int n_segments;
    std::vector<double> initial_segment_times;
    double initial_time_speed_ratio{0.65};  ///< warm-start segment time = seg_length / max(0.1, v_max * ratio)

    // ---- solver configuration ----
    int lbfgs_max_iterations{80};           ///< L-BFGS iteration cap

    // ---- motion camouflage ----
    bool use_vmc = true;                  ///< true=VMC降维(v参数), false=纯MINCO(直接优化xy)
    Eigen::Vector2d ref_point = Eigen::Vector2d::Zero();
    std::vector<Eigen::Vector2d> prey_points;
};

// ==================== Trajectory result ====================

/**
 * @struct McoTrajectory
 * @brief Output of the MCO optimiser
 */
struct McoTrajectory
{
    Eigen::MatrixXd coeffs;    ///< polynomial coefficients  (n_coeffs*n_segments) x dims
    Eigen::VectorXd T;         ///< per-segment durations     n_segments x 1
    double total_duration;     ///< sum of T
    double final_cost{0.0};    ///< cost at convergence

    bool isValid() const {
        return T.size() > 0 && coeffs.size() > 0;
    }
};

} // namespace mocha
