#pragma once

/**
 * @file prm_homotopy_planner.hpp
 * @brief Reference-guided visibility-PRM path planner using ESDF and path-level topology deduplication.
 *        Samples in a corridor around the current guide polyline, constructs a sparse
 *        guard/connector graph, enumerates simple paths, shortcuts them, smooths them
 *        while preserving topology, and deduplicates them by saturated H-signature.
 *        No ROS dependencies.
 */

#include <vector>
#include <string>
#include <functional>
#include <Eigen/Core>
#include "mocha_planner/core/types.hpp"

namespace mocha {

struct GeometricPath {
  std::vector<Eigen::Vector2d> points;
  std::string homotopy_signature;
  double length = 0.0;
};

struct PrmDebugStats {
  int samples_total = 0;
  int samples_accepted = 0;
  int guard_nodes = 0;
  int connector_nodes = 0;
  int start_degree = 0;
  int goal_degree = 0;
  int reachable_nodes_from_start = 0;
  bool goal_reachable_from_start = false;
  int graph_nodes = 0;
  int graph_edges = 0;
  int raw_simple_paths = 0;
  int unique_homotopy_paths = 0;
  int filtered_paths = 0;
};

class PrmHomotopyPlanner {
public:
  PrmHomotopyPlanner() = default;

  std::vector<GeometricPath> generatePaths(
      const Eigen::Vector2d& start,
      const Eigen::Vector2d& goal,
      const EsdfQueryFunc& esdf_query,
      const PrmParameters& prm_params,
      double robot_radius,
      const std::vector<Eigen::Vector2d>& obstacle_refs) const;

  void setPersistentSeedsFromPath(const std::vector<Eigen::Vector2d>& path);
  void clearPersistentSeeds();

  const PrmDebugStats& lastDebugStats() const { return last_debug_stats_; }

  static std::vector<Eigen::Vector2d> extractObstacleRefPoints(
      const std::vector<double>& sdf_grid,
      int width, int height,
      double resolution,
      double origin_x, double origin_y);

  static std::string computeHSignature(
      const std::vector<Eigen::Vector2d>& path,
      const std::vector<Eigen::Vector2d>& obstacle_refs);

  static std::vector<Eigen::Vector2d> resamplePolyline(
      const std::vector<Eigen::Vector2d>& path,
      double spacing);

  static bool edgeFree(
      const Eigen::Vector2d& a,
      const Eigen::Vector2d& b,
      const EsdfQueryFunc& esdf_query,
      double d_safe,
      int m_check);

private:
  mutable PrmDebugStats last_debug_stats_;
  mutable std::vector<Eigen::Vector2d> persistent_seed_nodes_;

  struct PrmGraph {
    std::vector<Eigen::Vector2d> nodes;
    std::vector<double> progress;
    std::vector<std::vector<std::pair<int, double>>> adj;
  };

  struct PathCandidate {
    std::vector<int> node_ids;
    double cost = 0.0;
  };

  /// Build a reference-guided visibility-PRM graph with sparse guard/connector nodes.
  static PrmGraph buildPrm(
      const Eigen::Vector2d& start,
      const Eigen::Vector2d& goal,
      const EsdfQueryFunc& esdf_query,
      const PrmParameters& params,
      double d_safe,
      const std::vector<Eigen::Vector2d>& propagated_samples,
      const std::vector<Eigen::Vector2d>& relevant_obstacle_refs,
      PrmDebugStats* debug_stats);

  static double pathLength(
      const PrmGraph& graph,
      const std::vector<int>& node_ids);
  static bool hasEdge(
      const PrmGraph& graph,
      int u,
      int v);
  static std::vector<PathCandidate> enumerateLooplessPaths(
      const PrmGraph& graph,
      int src,
      int dst,
      int max_candidates,
      int max_depth);
};

} // namespace mocha
