#include "mocha_planner/core/prm_homotopy_planner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <random>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mocha {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kCollisionCheckStep = 0.10;
constexpr double kGeometryEps = 1e-9;
constexpr double kProgressEps = 1e-4;
constexpr int kMaxSearchMultiplier = 128;

double clampAngle(double angle)
{
  while (angle <= -kPi) {
    angle += 2.0 * kPi;
  }
  while (angle > kPi) {
    angle -= 2.0 * kPi;
  }
  return angle;
}

double cross2d(const Eigen::Vector2d& a, const Eigen::Vector2d& b)
{
  return a.x() * b.y() - a.y() * b.x();
}

double pointToSegmentDistance(
    const Eigen::Vector2d& point,
    const Eigen::Vector2d& seg_a,
    const Eigen::Vector2d& seg_b)
{
  const Eigen::Vector2d ab = seg_b - seg_a;
  const double denom = ab.squaredNorm();
  if (denom < kGeometryEps) {
    return (point - seg_a).norm();
  }

  const double t = std::clamp((point - seg_a).dot(ab) / denom, 0.0, 1.0);
  const Eigen::Vector2d projection = seg_a + t * ab;
  return (point - projection).norm();
}

double pointToPolylineDistance(
    const Eigen::Vector2d& point,
    const std::vector<Eigen::Vector2d>& polyline)
{
  if (polyline.empty()) {
    return std::numeric_limits<double>::infinity();
  }
  if (polyline.size() == 1) {
    return (point - polyline.front()).norm();
  }

  double best = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < polyline.size(); ++i) {
    best = std::min(best, pointToSegmentDistance(point, polyline[i - 1], polyline[i]));
  }
  return best;
}

double polylineLength(const std::vector<Eigen::Vector2d>& polyline)
{
  double length = 0.0;
  for (size_t i = 1; i < polyline.size(); ++i) {
    length += (polyline[i] - polyline[i - 1]).norm();
  }
  return length;
}

double projectPointToPolylineProgress(
    const Eigen::Vector2d& point,
    const std::vector<Eigen::Vector2d>& polyline)
{
  if (polyline.empty()) {
    return 0.0;
  }
  if (polyline.size() == 1) {
    return 0.0;
  }

  double accumulated_length = 0.0;
  double best_progress = 0.0;
  double best_distance_sq = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < polyline.size(); ++i) {
    const Eigen::Vector2d segment = polyline[i] - polyline[i - 1];
    const double segment_length_sq = segment.squaredNorm();
    const double segment_length = std::sqrt(segment_length_sq);
    double alpha = 0.0;
    if (segment_length_sq > kGeometryEps) {
      alpha = std::clamp(
          (point - polyline[i - 1]).dot(segment) / segment_length_sq,
          0.0,
          1.0);
    }
    const Eigen::Vector2d projection =
        polyline[i - 1] + alpha * segment;
    const double distance_sq = (point - projection).squaredNorm();
    if (distance_sq < best_distance_sq) {
      best_distance_sq = distance_sq;
      best_progress = accumulated_length + alpha * segment_length;
    }
    accumulated_length += segment_length;
  }

  return best_progress;
}

bool isDuplicateNode(
    const std::vector<Eigen::Vector2d>& nodes,
    const Eigen::Vector2d& candidate,
    double min_spacing)
{
  const double min_spacing_sq = min_spacing * min_spacing;
  for (const auto& node : nodes) {
    if ((node - candidate).squaredNorm() <= min_spacing_sq) {
      return true;
    }
  }
  return false;
}

int orientation(
    const Eigen::Vector2d& a,
    const Eigen::Vector2d& b,
    const Eigen::Vector2d& c)
{
  const double value = cross2d(b - a, c - a);
  if (value > kGeometryEps) {
    return 1;
  }
  if (value < -kGeometryEps) {
    return -1;
  }
  return 0;
}

bool onSegment(
    const Eigen::Vector2d& a,
    const Eigen::Vector2d& b,
    const Eigen::Vector2d& p)
{
  return
      p.x() >= std::min(a.x(), b.x()) - kGeometryEps &&
      p.x() <= std::max(a.x(), b.x()) + kGeometryEps &&
      p.y() >= std::min(a.y(), b.y()) - kGeometryEps &&
      p.y() <= std::max(a.y(), b.y()) + kGeometryEps;
}

bool segmentsIntersect(
    const Eigen::Vector2d& a,
    const Eigen::Vector2d& b,
    const Eigen::Vector2d& c,
    const Eigen::Vector2d& d)
{
  const int o1 = orientation(a, b, c);
  const int o2 = orientation(a, b, d);
  const int o3 = orientation(c, d, a);
  const int o4 = orientation(c, d, b);

  if (o1 != o2 && o3 != o4) {
    return true;
  }
  if (o1 == 0 && onSegment(a, b, c)) {
    return true;
  }
  if (o2 == 0 && onSegment(a, b, d)) {
    return true;
  }
  if (o3 == 0 && onSegment(c, d, a)) {
    return true;
  }
  if (o4 == 0 && onSegment(c, d, b)) {
    return true;
  }
  return false;
}

bool hasSelfIntersection(const std::vector<Eigen::Vector2d>& polyline)
{
  if (polyline.size() < 4) {
    return false;
  }

  for (size_t i = 1; i < polyline.size(); ++i) {
    const Eigen::Vector2d& a0 = polyline[i - 1];
    const Eigen::Vector2d& a1 = polyline[i];
    for (size_t j = i + 2; j < polyline.size(); ++j) {
      const Eigen::Vector2d& b0 = polyline[j - 1];
      const Eigen::Vector2d& b1 = polyline[j];
      if (segmentsIntersect(a0, a1, b0, b1)) {
        return true;
      }
    }
  }

  return false;
}

std::vector<int> computeRawSignatureVector(
    const std::vector<Eigen::Vector2d>& path,
    const std::vector<Eigen::Vector2d>& obstacle_refs)
{
  std::vector<int> signature(obstacle_refs.size(), 0);
  if (path.size() < 2) {
    return signature;
  }

  for (size_t k = 0; k < obstacle_refs.size(); ++k) {
    const Eigen::Vector2d& obs = obstacle_refs[k];
    const Eigen::Vector2d start_vec = path.front() - obs;
    const Eigen::Vector2d goal_vec = path.back() - obs;
    const double endpoint_delta = clampAngle(
        std::atan2(goal_vec.y(), goal_vec.x()) -
        std::atan2(start_vec.y(), start_vec.x()));

    double winding = 0.0;
    for (size_t i = 1; i < path.size(); ++i) {
      const Eigen::Vector2d a = path[i - 1] - obs;
      const Eigen::Vector2d b = path[i] - obs;
      const double angle_a = std::atan2(a.y(), a.x());
      const double angle_b = std::atan2(b.y(), b.x());
      winding += clampAngle(angle_b - angle_a);
    }

    signature[k] = static_cast<int>(
        std::llround((winding - endpoint_delta) / (2.0 * kPi)));
  }

  return signature;
}

std::vector<int> saturateSignature(const std::vector<int>& raw_signature)
{
  std::vector<int> saturated(raw_signature.size(), 0);
  for (size_t i = 0; i < raw_signature.size(); ++i) {
    saturated[i] = (raw_signature[i] > 0) - (raw_signature[i] < 0);
  }
  return saturated;
}

std::string signatureKey(const std::vector<int>& signature)
{
  std::ostringstream oss;
  bool first = true;
  for (size_t i = 0; i < signature.size(); ++i) {
    if (signature[i] == 0) {
      continue;
    }
    if (!first) {
      oss << ";";
    }
    first = false;
    oss << "b[" << i << "]=" << signature[i];
  }

  const std::string key = oss.str();
  return key.empty() ? "w=0" : key;
}

double lateralSamplingWidth(
    const PrmParameters& params,
    double d_safe)
{
  const double lower = std::max(2.0 * d_safe, 0.25);
  const double upper = std::max(lower, params.local_range);
  const double nominal = params.lateral_sample_factor * upper;
  return std::clamp(nominal, lower, upper);
}

struct PolylineSampler
{
  explicit PolylineSampler(std::vector<Eigen::Vector2d> polyline)
  : polyline_(std::move(polyline))
  {
    if (polyline_.empty()) {
      return;
    }
    cumulative_length_.resize(polyline_.size(), 0.0);
    for (size_t i = 1; i < polyline_.size(); ++i) {
      cumulative_length_[i] =
          cumulative_length_[i - 1] + (polyline_[i] - polyline_[i - 1]).norm();
    }
    total_length_ = cumulative_length_.back();
  }

  Eigen::Vector2d sample(double arc_length, Eigen::Vector2d* tangent) const
  {
    if (polyline_.empty()) {
      if (tangent != nullptr) {
        *tangent = Eigen::Vector2d::UnitX();
      }
      return Eigen::Vector2d::Zero();
    }
    if (polyline_.size() == 1 || total_length_ < kGeometryEps) {
      if (tangent != nullptr) {
        *tangent = Eigen::Vector2d::UnitX();
      }
      return polyline_.front();
    }

    const double clamped_s = std::clamp(arc_length, 0.0, total_length_);
    auto upper = std::lower_bound(cumulative_length_.begin(), cumulative_length_.end(), clamped_s);
    size_t seg_index = static_cast<size_t>(std::distance(cumulative_length_.begin(), upper));
    if (seg_index == 0) {
      seg_index = 1;
    }
    if (seg_index >= polyline_.size()) {
      seg_index = polyline_.size() - 1;
    }

    const double seg_start_s = cumulative_length_[seg_index - 1];
    const double seg_end_s = cumulative_length_[seg_index];
    const double seg_length = std::max(seg_end_s - seg_start_s, kGeometryEps);
    const double alpha = (clamped_s - seg_start_s) / seg_length;
    const Eigen::Vector2d direction = (polyline_[seg_index] - polyline_[seg_index - 1]).normalized();
    if (tangent != nullptr) {
      *tangent = direction;
    }
    return (1.0 - alpha) * polyline_[seg_index - 1] + alpha * polyline_[seg_index];
  }

  double totalLength() const
  {
    return total_length_;
  }

private:
  std::vector<Eigen::Vector2d> polyline_;
  std::vector<double> cumulative_length_;
  double total_length_{0.0};
};

std::vector<Eigen::Vector2d> sanitizeGuidePolyline(
    const Eigen::Vector2d& start,
    const Eigen::Vector2d& goal,
    const std::vector<Eigen::Vector2d>& guide_points)
{
  std::vector<Eigen::Vector2d> guide = guide_points;
  if (guide.empty()) {
    guide = {start, goal};
  }

  if ((guide.front() - start).norm() > 1e-6) {
    guide.insert(guide.begin(), start);
  } else {
    guide.front() = start;
  }
  if ((guide.back() - goal).norm() > 1e-6) {
    guide.push_back(goal);
  } else {
    guide.back() = goal;
  }

  std::vector<Eigen::Vector2d> compacted;
  compacted.reserve(guide.size());
  for (const auto& point : guide) {
    if (compacted.empty() || (compacted.back() - point).norm() > 1e-6) {
      compacted.push_back(point);
    }
  }
  if (compacted.size() == 1) {
    compacted.push_back(goal);
  }
  return compacted;
}

std::vector<Eigen::Vector2d> selectRelevantObstacleRefs(
    const std::vector<Eigen::Vector2d>& guide_polyline,
    const std::vector<Eigen::Vector2d>& obstacle_refs,
    double corridor_margin,
    size_t max_refs)
{
  if (obstacle_refs.empty()) {
    return {};
  }

  std::vector<std::pair<double, Eigen::Vector2d>> ranked_refs;
  ranked_refs.reserve(obstacle_refs.size());
  for (const auto& ref : obstacle_refs) {
    ranked_refs.emplace_back(pointToPolylineDistance(ref, guide_polyline), ref);
  }

  std::sort(
      ranked_refs.begin(),
      ranked_refs.end(),
      [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

  std::vector<Eigen::Vector2d> selected;
  selected.reserve(std::min(max_refs, ranked_refs.size()));
  for (const auto& [distance, ref] : ranked_refs) {
    if (distance <= corridor_margin || selected.size() < max_refs) {
      selected.push_back(ref);
    }
    if (selected.size() >= max_refs) {
      break;
    }
  }
  return selected;
}

std::uint64_t anchorPairKey(int a, int b)
{
  if (a > b) {
    std::swap(a, b);
  }
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
         static_cast<std::uint32_t>(b);
}

bool sampleProgressStrictlyBetween(
    double sample_progress,
    double progress_a,
    double progress_b)
{
  const double lower = std::min(progress_a, progress_b);
  const double upper = std::max(progress_a, progress_b);
  return
      sample_progress > lower + kProgressEps &&
      sample_progress < upper - kProgressEps;
}

bool projectPointToSafeClearance(
    Eigen::Vector2d* point,
    const EsdfQueryFunc& esdf_query,
    double d_safe)
{
  double distance = 0.0;
  Eigen::Vector2d gradient = Eigen::Vector2d::Zero();
  if (!esdf_query(*point, distance, gradient)) {
    return false;
  }
  if (distance >= d_safe) {
    return true;
  }

  for (int i = 0; i < 24 && distance < d_safe; ++i) {
    if (gradient.norm() < 1e-6) {
      break;
    }
    *point += gradient.normalized() * std::max(0.02, d_safe - distance);
    if (!esdf_query(*point, distance, gradient)) {
      return false;
    }
  }

  return distance >= d_safe;
}

std::vector<Eigen::Vector2d> shortcutHomotopyPath(
    const std::vector<Eigen::Vector2d>& path,
    const std::vector<Eigen::Vector2d>& obstacle_refs,
    const EsdfQueryFunc& esdf_query,
    double d_safe,
    int m_check)
{
  if (path.size() < 3) {
    return path;
  }

  const std::vector<int> original_signature =
      computeRawSignatureVector(path, obstacle_refs);
  std::vector<Eigen::Vector2d> shortened = path;

  bool changed = true;
  while (changed && shortened.size() >= 3) {
    changed = false;
    for (size_t i = 0; i + 2 < shortened.size() && !changed; ++i) {
      for (size_t j = shortened.size() - 1; j > i + 1; --j) {
        if (!PrmHomotopyPlanner::edgeFree(
                shortened[i],
                shortened[j],
                esdf_query,
                d_safe,
                m_check)) {
          continue;
        }

        std::vector<Eigen::Vector2d> candidate;
        candidate.reserve(shortened.size() - (j - i - 1));
        candidate.insert(candidate.end(), shortened.begin(), shortened.begin() + static_cast<long>(i + 1));
        candidate.insert(candidate.end(), shortened.begin() + static_cast<long>(j), shortened.end());

        if (hasSelfIntersection(candidate)) {
          continue;
        }
        if (computeRawSignatureVector(candidate, obstacle_refs) != original_signature) {
          continue;
        }

        shortened = std::move(candidate);
        changed = true;
        break;
      }
    }
  }

  return shortened;
}

bool polylineCollisionFree(
    const std::vector<Eigen::Vector2d>& path,
    const EsdfQueryFunc& esdf_query,
    double d_safe,
    int m_check)
{
  if (path.size() < 2) {
    return false;
  }

  for (size_t i = 1; i < path.size(); ++i) {
    if (!PrmHomotopyPlanner::edgeFree(
            path[i - 1],
            path[i],
            esdf_query,
            d_safe,
            m_check)) {
      return false;
    }
  }
  return true;
}

std::vector<Eigen::Vector2d> smoothHomotopyPath(
    const std::vector<Eigen::Vector2d>& path,
    const std::vector<Eigen::Vector2d>& obstacle_refs,
    const EsdfQueryFunc& esdf_query,
    double d_safe,
    int m_check,
    double smooth_spacing,
    int smooth_passes,
    double smooth_strength)
{
  if (path.size() < 3 || smooth_passes <= 0 || smooth_strength <= 0.0) {
    return path;
  }

  std::vector<Eigen::Vector2d> smoothed =
      PrmHomotopyPlanner::resamplePolyline(path, std::max(0.05, smooth_spacing));
  if (smoothed.size() < 3) {
    return path;
  }

  const std::vector<int> original_signature = computeRawSignatureVector(smoothed, obstacle_refs);
  smoothed.front() = path.front();
  smoothed.back() = path.back();

  for (int pass = 0; pass < smooth_passes; ++pass) {
    std::vector<Eigen::Vector2d> candidate = smoothed;
    bool moved = false;

    for (size_t i = 1; i + 1 < smoothed.size(); ++i) {
      const Eigen::Vector2d laplacian =
          0.5 * (smoothed[i - 1] + smoothed[i + 1]) - smoothed[i];
      if (laplacian.squaredNorm() < 1e-10) {
        continue;
      }

      Eigen::Vector2d proposal = smoothed[i] + smooth_strength * laplacian;
      if (!projectPointToSafeClearance(&proposal, esdf_query, d_safe)) {
        continue;
      }

      candidate[i] = proposal;
      moved = true;
    }

    if (!moved) {
      break;
    }

    candidate.front() = path.front();
    candidate.back() = path.back();
    if (hasSelfIntersection(candidate) ||
        !polylineCollisionFree(candidate, esdf_query, d_safe, m_check) ||
        computeRawSignatureVector(candidate, obstacle_refs) != original_signature) {
      continue;
    }

    smoothed = std::move(candidate);
  }

  return smoothed;
}

}  // namespace

void PrmHomotopyPlanner::setPersistentSeedsFromPath(const std::vector<Eigen::Vector2d>& path)
{
  persistent_seed_nodes_.clear();
  if (path.size() <= 2) {
    return;
  }

  persistent_seed_nodes_.reserve(path.size() - 2);
  for (size_t i = 1; i + 1 < path.size(); ++i) {
    persistent_seed_nodes_.push_back(path[i]);
  }
}

void PrmHomotopyPlanner::clearPersistentSeeds()
{
  persistent_seed_nodes_.clear();
}

std::vector<Eigen::Vector2d> PrmHomotopyPlanner::extractObstacleRefPoints(
    const std::vector<double>& sdf_grid,
    int width,
    int height,
    double resolution,
    double origin_x,
    double origin_y)
{
  std::vector<Eigen::Vector2d> refs;
  const int cell_count = width * height;
  if (static_cast<int>(sdf_grid.size()) < cell_count) {
    return refs;
  }

  std::vector<bool> visited(static_cast<size_t>(cell_count), false);
  const int dx8[] = {-1, 0, 1, -1, 1, -1, 0, 1};
  const int dy8[] = {-1, -1, -1, 0, 0, 1, 1, 1};

  for (int idx = 0; idx < cell_count; ++idx) {
    if (visited[static_cast<size_t>(idx)] || sdf_grid[static_cast<size_t>(idx)] >= 0.0) {
      continue;
    }

    std::queue<int> bfs;
    bfs.push(idx);
    visited[static_cast<size_t>(idx)] = true;

    double sum_x = 0.0;
    double sum_y = 0.0;
    int count = 0;

    while (!bfs.empty()) {
      const int current = bfs.front();
      bfs.pop();
      const int cx = current % width;
      const int cy = current / width;

      sum_x += origin_x + (static_cast<double>(cx) + 0.5) * resolution;
      sum_y += origin_y + (static_cast<double>(cy) + 0.5) * resolution;
      ++count;

      for (int d = 0; d < 8; ++d) {
        const int nx = cx + dx8[d];
        const int ny = cy + dy8[d];
        if (nx < 0 || nx >= width || ny < 0 || ny >= height) {
          continue;
        }
        const int next = ny * width + nx;
        if (!visited[static_cast<size_t>(next)] && sdf_grid[static_cast<size_t>(next)] < 0.0) {
          visited[static_cast<size_t>(next)] = true;
          bfs.push(next);
        }
      }
    }

    if (count > 0) {
      refs.emplace_back(sum_x / count, sum_y / count);
    }
  }

  return refs;
}

std::string PrmHomotopyPlanner::computeHSignature(
    const std::vector<Eigen::Vector2d>& path,
    const std::vector<Eigen::Vector2d>& obstacle_refs)
{
  return signatureKey(computeRawSignatureVector(path, obstacle_refs));
}

bool PrmHomotopyPlanner::edgeFree(
    const Eigen::Vector2d& a,
    const Eigen::Vector2d& b,
    const EsdfQueryFunc& esdf_query,
    double d_safe,
    int m_check)
{
  const double edge_length = (b - a).norm();
  const int sampling_steps = std::max(
      std::max(1, m_check),
      static_cast<int>(std::ceil(edge_length / kCollisionCheckStep)));

  for (int step = 0; step <= sampling_steps; ++step) {
    const double t = static_cast<double>(step) / static_cast<double>(sampling_steps);
    const Eigen::Vector2d point = (1.0 - t) * a + t * b;
    double distance = 0.0;
    Eigen::Vector2d gradient = Eigen::Vector2d::Zero();
    if (!esdf_query(point, distance, gradient) || distance < d_safe) {
      return false;
    }
  }

  return true;
}

std::vector<Eigen::Vector2d> PrmHomotopyPlanner::resamplePolyline(
    const std::vector<Eigen::Vector2d>& path,
    double spacing)
{
  if (path.size() < 2 || spacing <= 0.0) {
    return path;
  }

  std::vector<double> cumulative_length(path.size(), 0.0);
  for (size_t i = 1; i < path.size(); ++i) {
    cumulative_length[i] = cumulative_length[i - 1] + (path[i] - path[i - 1]).norm();
  }

  const double total_length = cumulative_length.back();
  if (total_length < kGeometryEps) {
    return {path.front()};
  }

  const int point_count = std::max(2, static_cast<int>(std::ceil(total_length / spacing)) + 1);
  const double actual_spacing = total_length / static_cast<double>(point_count - 1);

  std::vector<Eigen::Vector2d> result;
  result.reserve(static_cast<size_t>(point_count));
  result.push_back(path.front());

  size_t segment_index = 0;
  for (int i = 1; i < point_count - 1; ++i) {
    const double target_length = actual_spacing * static_cast<double>(i);
    while (segment_index + 1 < path.size() - 1 &&
           cumulative_length[segment_index + 1] < target_length) {
      ++segment_index;
    }
    const double segment_length =
        cumulative_length[segment_index + 1] - cumulative_length[segment_index];
    const double alpha = (segment_length > kGeometryEps)
        ? (target_length - cumulative_length[segment_index]) / segment_length
        : 0.0;
    result.push_back(
        (1.0 - alpha) * path[segment_index] + alpha * path[segment_index + 1]);
  }

  result.push_back(path.back());
  return result;
}

double PrmHomotopyPlanner::pathLength(
    const PrmGraph& graph,
    const std::vector<int>& node_ids)
{
  double length = 0.0;
  for (size_t i = 1; i < node_ids.size(); ++i) {
    length += (graph.nodes[node_ids[i]] - graph.nodes[node_ids[i - 1]]).norm();
  }
  return length;
}

bool PrmHomotopyPlanner::hasEdge(
    const PrmGraph& graph,
    int u,
    int v)
{
  for (const auto& [next, _weight] : graph.adj[static_cast<size_t>(u)]) {
    if (next == v) {
      return true;
    }
  }
  return false;
}

PrmHomotopyPlanner::PrmGraph PrmHomotopyPlanner::buildPrm(
    const Eigen::Vector2d& start,
    const Eigen::Vector2d& goal,
    const EsdfQueryFunc& esdf_query,
    const PrmParameters& params,
    double d_safe,
    const std::vector<Eigen::Vector2d>& propagated_samples,
    const std::vector<Eigen::Vector2d>& relevant_obstacle_refs,
    PrmDebugStats* debug_stats)
{
  const std::vector<Eigen::Vector2d> guide_polyline = sanitizeGuidePolyline(start, goal, params.guide_points);
  const PolylineSampler guide_sampler(guide_polyline);
  const double total_guide_length = std::max(guide_sampler.totalLength(), kGeometryEps);
  const double sample_width = lateralSamplingWidth(params, d_safe);

  PrmGraph graph;
  graph.nodes = {start, goal};
  graph.progress = {0.0, total_guide_length};
  graph.adj.resize(2);

  std::vector<bool> is_guard = {true, false};
  const double min_guard_spacing =
      params.min_guard_spacing > 0.0 ? params.min_guard_spacing : std::max(2.0 * d_safe, 0.35);
  const double min_sample_spacing = std::max(0.5 * d_safe, 0.10);

  auto addUndirectedEdge = [&graph](int u, int v) {
    if (u == v || PrmHomotopyPlanner::hasEdge(graph, u, v)) {
      return;
    }
    const double weight = (graph.nodes[static_cast<size_t>(u)] - graph.nodes[static_cast<size_t>(v)]).norm();
    graph.adj[static_cast<size_t>(u)].push_back({v, weight});
    graph.adj[static_cast<size_t>(v)].push_back({u, weight});
  };

  auto updateUndirectedEdge = [&graph](int u, int v) {
    const double weight = (graph.nodes[static_cast<size_t>(u)] - graph.nodes[static_cast<size_t>(v)]).norm();
    for (auto& [node_id, edge_weight] : graph.adj[static_cast<size_t>(u)]) {
      if (node_id == v) {
        edge_weight = weight;
        break;
      }
    }
    for (auto& [node_id, edge_weight] : graph.adj[static_cast<size_t>(v)]) {
      if (node_id == u) {
        edge_weight = weight;
        break;
      }
    }
  };

  auto findVisibleGuards = [&](const Eigen::Vector2d& sample) {
    std::vector<int> visible_guards;
    for (size_t node_id = 0; node_id < graph.nodes.size(); ++node_id) {
      if (!is_guard[node_id]) {
        continue;
      }
      if (edgeFree(sample, graph.nodes[node_id], esdf_query, d_safe, params.m_check)) {
        visible_guards.push_back(static_cast<int>(node_id));
      }
    }
    return visible_guards;
  };

  std::unordered_map<std::uint64_t, std::unordered_map<std::string, int>> connector_lookup;
  auto registerConnector = [&](int anchor_a, int anchor_b, const Eigen::Vector2d& sample, double sample_progress) {
    if (anchor_a == anchor_b) {
      return;
    }
    if (!sampleProgressStrictlyBetween(
            sample_progress,
            graph.progress[static_cast<size_t>(anchor_a)],
            graph.progress[static_cast<size_t>(anchor_b)])) {
      return;
    }

    int forward_anchor = anchor_a;
    int backward_anchor = anchor_b;
    if (graph.progress[static_cast<size_t>(forward_anchor)] >
        graph.progress[static_cast<size_t>(backward_anchor)]) {
      std::swap(forward_anchor, backward_anchor);
    }

    const std::vector<Eigen::Vector2d> segment = {
      graph.nodes[static_cast<size_t>(forward_anchor)],
      sample,
      graph.nodes[static_cast<size_t>(backward_anchor)]
    };
    const std::string topology_key = signatureKey(
        saturateSignature(computeRawSignatureVector(segment, relevant_obstacle_refs)));
    const std::uint64_t pair_key = anchorPairKey(anchor_a, anchor_b);
    const double connector_cost =
        (graph.nodes[static_cast<size_t>(forward_anchor)] - sample).norm() +
        (sample - graph.nodes[static_cast<size_t>(backward_anchor)]).norm();

    auto& pair_entries = connector_lookup[pair_key];
    const auto existing = pair_entries.find(topology_key);
    if (existing != pair_entries.end()) {
      const int connector_id = existing->second;
      const double old_cost =
          (graph.nodes[static_cast<size_t>(forward_anchor)] - graph.nodes[static_cast<size_t>(connector_id)]).norm() +
          (graph.nodes[static_cast<size_t>(connector_id)] - graph.nodes[static_cast<size_t>(backward_anchor)]).norm();
      if (connector_cost + 1e-6 < old_cost) {
        graph.nodes[static_cast<size_t>(connector_id)] = sample;
        graph.progress[static_cast<size_t>(connector_id)] = sample_progress;
        updateUndirectedEdge(connector_id, anchor_a);
        updateUndirectedEdge(connector_id, anchor_b);
      }
      return;
    }

    const int connector_id = static_cast<int>(graph.nodes.size());
    graph.nodes.push_back(sample);
    graph.progress.push_back(sample_progress);
    graph.adj.emplace_back();
    is_guard.push_back(false);
    addUndirectedEdge(connector_id, anchor_a);
    addUndirectedEdge(connector_id, anchor_b);
    pair_entries[topology_key] = connector_id;
  };

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> arc_dist(
      0.0,
      total_guide_length);
  std::uniform_real_distribution<double> lateral_dist(-sample_width, sample_width);

  auto tryInsertSample = [&](Eigen::Vector2d sample, bool boundary_bias) {
    double distance = 0.0;
    Eigen::Vector2d gradient = Eigen::Vector2d::Zero();
    if (!esdf_query(sample, distance, gradient)) {
      return;
    }

    if (boundary_bias) {
      const double target_distance = 1.25 * d_safe;
      if (distance > target_distance && gradient.norm() > 1e-6) {
        sample -= gradient.normalized() * (distance - target_distance);
        if (!esdf_query(sample, distance, gradient)) {
          return;
        }
      }
    }

    if (distance < d_safe && !projectPointToSafeClearance(&sample, esdf_query, d_safe)) {
      return;
    }

    if (debug_stats != nullptr) {
      ++debug_stats->samples_accepted;
    }

    const double sample_progress = projectPointToPolylineProgress(sample, guide_polyline);
    if (sample_progress <= kProgressEps || sample_progress >= total_guide_length - kProgressEps) {
      return;
    }
    if (isDuplicateNode(graph.nodes, sample, min_sample_spacing)) {
      return;
    }

    const bool goal_visible = edgeFree(sample, goal, esdf_query, d_safe, params.m_check);
    const std::vector<int> visible_guards = findVisibleGuards(sample);

    if (!goal_visible && visible_guards.empty()) {
      std::vector<Eigen::Vector2d> guard_positions;
      guard_positions.reserve(graph.nodes.size());
      for (size_t node_id = 0; node_id < graph.nodes.size(); ++node_id) {
        if (is_guard[node_id]) {
          guard_positions.push_back(graph.nodes[node_id]);
        }
      }
      if (!isDuplicateNode(guard_positions, sample, min_guard_spacing)) {
        graph.nodes.push_back(sample);
        graph.progress.push_back(sample_progress);
        graph.adj.emplace_back();
        is_guard.push_back(true);
      }
      return;
    }

    if (!goal_visible && visible_guards.size() == 2) {
      registerConnector(visible_guards[0], visible_guards[1], sample, sample_progress);
      return;
    }

    if (goal_visible && visible_guards.size() == 1) {
      registerConnector(visible_guards[0], 1, sample, sample_progress);
    }
  };

  for (const auto& seed : propagated_samples) {
    tryInsertSample(seed, false);
  }

  for (int sample_index = 0; sample_index < params.n_samples; ++sample_index) {
    if (debug_stats != nullptr) {
      ++debug_stats->samples_total;
    }

    Eigen::Vector2d tangent = Eigen::Vector2d::UnitX();
    const double guide_progress = arc_dist(rng);
    const Eigen::Vector2d guide_center = guide_sampler.sample(guide_progress, &tangent);
    if (tangent.squaredNorm() < kGeometryEps) {
      tangent = Eigen::Vector2d::UnitX();
    }
    const Eigen::Vector2d normal(-tangent.y(), tangent.x());
    Eigen::Vector2d sample = guide_center + lateral_dist(rng) * normal;
    tryInsertSample(
        sample,
        sample_index < static_cast<int>(params.n_samples * params.boundary_sample_ratio));
  }

  if (debug_stats != nullptr) {
    debug_stats->graph_nodes = static_cast<int>(graph.nodes.size());
    debug_stats->guard_nodes = 0;
    debug_stats->connector_nodes = 0;
    for (size_t node_id = 0; node_id < graph.nodes.size(); ++node_id) {
      if (node_id == 1) {
        continue;
      }
      if (is_guard[node_id]) {
        ++debug_stats->guard_nodes;
      } else {
        ++debug_stats->connector_nodes;
      }
    }

    int total_edges = 0;
    for (const auto& neighbors : graph.adj) {
      total_edges += static_cast<int>(neighbors.size());
    }
    debug_stats->graph_edges = total_edges / 2;
    debug_stats->start_degree = static_cast<int>(graph.adj.front().size());
    debug_stats->goal_degree = static_cast<int>(graph.adj[1].size());

    std::vector<bool> visited(graph.nodes.size(), false);
    std::queue<int> bfs;
    bfs.push(0);
    visited[0] = true;
    while (!bfs.empty()) {
      const int current = bfs.front();
      bfs.pop();
      ++debug_stats->reachable_nodes_from_start;
      for (const auto& [next, _weight] : graph.adj[static_cast<size_t>(current)]) {
        if (!visited[static_cast<size_t>(next)]) {
          visited[static_cast<size_t>(next)] = true;
          bfs.push(next);
        }
      }
    }
    debug_stats->goal_reachable_from_start =
        visited.size() > 1 ? visited[1] : false;
  }

  return graph;
}

std::vector<PrmHomotopyPlanner::PathCandidate> PrmHomotopyPlanner::enumerateLooplessPaths(
    const PrmGraph& graph,
    int src,
    int dst,
    int max_candidates,
    int max_depth)
{
  std::vector<PathCandidate> results;
  if (max_candidates <= 0 || src < 0 || dst < 0 || src >= static_cast<int>(graph.nodes.size()) ||
      dst >= static_cast<int>(graph.nodes.size())) {
    return results;
  }

  struct PartialPath
  {
    std::vector<int> node_ids;
    double cost{0.0};
    double score{0.0};
  };

  struct Compare
  {
    bool operator()(const PartialPath& lhs, const PartialPath& rhs) const
    {
      if (lhs.score == rhs.score) {
        return lhs.cost > rhs.cost;
      }
      return lhs.score > rhs.score;
    }
  };

  const auto heuristic = [&graph, dst](int node_id) {
    return (graph.nodes[static_cast<size_t>(node_id)] - graph.nodes[static_cast<size_t>(dst)]).norm();
  };
  const auto advancesAlongGuide = [&graph](int from, int to) {
    return
        graph.progress[static_cast<size_t>(to)] >
        graph.progress[static_cast<size_t>(from)] + kProgressEps;
  };

  std::priority_queue<PartialPath, std::vector<PartialPath>, Compare> queue;
  queue.push(PartialPath{{src}, 0.0, heuristic(src)});

  int expansions = 0;
  const int max_expansions = std::max(max_candidates * kMaxSearchMultiplier, 256);

  while (!queue.empty() && static_cast<int>(results.size()) < max_candidates && expansions < max_expansions) {
    PartialPath current = queue.top();
    queue.pop();
    ++expansions;

    const int node_id = current.node_ids.back();
    if (node_id == dst) {
      results.push_back(PathCandidate{current.node_ids, current.cost});
      continue;
    }

    if (static_cast<int>(current.node_ids.size()) >= max_depth) {
      continue;
    }

    std::vector<std::pair<int, double>> neighbors = graph.adj[static_cast<size_t>(node_id)];
    std::sort(
        neighbors.begin(),
        neighbors.end(),
        [&heuristic](const auto& lhs, const auto& rhs) {
          return lhs.second + heuristic(lhs.first) < rhs.second + heuristic(rhs.first);
        });

    for (const auto& [next, weight] : neighbors) {
      if (std::find(current.node_ids.begin(), current.node_ids.end(), next) != current.node_ids.end()) {
        continue;
      }
      if (!advancesAlongGuide(node_id, next)) {
        continue;
      }

      PartialPath extended = current;
      extended.node_ids.push_back(next);
      extended.cost += weight;
      extended.score = extended.cost + heuristic(next);
      queue.push(std::move(extended));
    }
  }

  std::sort(
      results.begin(),
      results.end(),
      [](const PathCandidate& lhs, const PathCandidate& rhs) { return lhs.cost < rhs.cost; });
  return results;
}

std::vector<GeometricPath> PrmHomotopyPlanner::generatePaths(
    const Eigen::Vector2d& start,
    const Eigen::Vector2d& goal,
    const EsdfQueryFunc& esdf_query,
    const PrmParameters& prm_params,
    double robot_radius,
    const std::vector<Eigen::Vector2d>& obstacle_refs) const
{
  last_debug_stats_ = PrmDebugStats{};

  const double d_safe = prm_params.d_safe > 0.0 ? prm_params.d_safe : robot_radius + 0.1;
  const std::vector<Eigen::Vector2d> guide_polyline = sanitizeGuidePolyline(start, goal, prm_params.guide_points);
  const double corridor_margin =
      lateralSamplingWidth(prm_params, d_safe) + 2.0 * d_safe;
  const std::vector<Eigen::Vector2d> relevant_obstacle_refs = selectRelevantObstacleRefs(
      guide_polyline,
      obstacle_refs,
      corridor_margin,
      std::max(1, prm_params.max_obstacle_refs));

  PrmGraph graph = buildPrm(
      start,
      goal,
      esdf_query,
      prm_params,
      d_safe,
      persistent_seed_nodes_,
      relevant_obstacle_refs,
      &last_debug_stats_);
  if (graph.nodes.size() < 2) {
    return {};
  }

  const std::vector<PathCandidate> raw_paths = enumerateLooplessPaths(
      graph,
      0,
      1,
      std::max(prm_params.k_max_paths, prm_params.max_raw_paths),
      std::max(2, prm_params.max_path_depth));
  last_debug_stats_.raw_simple_paths = static_cast<int>(raw_paths.size());

  std::unordered_map<std::string, GeometricPath> best_by_signature;
  for (const auto& candidate : raw_paths) {
    GeometricPath path;
    path.points.reserve(candidate.node_ids.size());
    for (int node_id : candidate.node_ids) {
      path.points.push_back(graph.nodes[static_cast<size_t>(node_id)]);
    }
    if (path.points.size() < 2 || hasSelfIntersection(path.points)) {
      continue;
    }

    path.points = shortcutHomotopyPath(
        path.points,
        relevant_obstacle_refs,
        esdf_query,
        d_safe,
        prm_params.m_check);
    if (path.points.size() < 2 || hasSelfIntersection(path.points)) {
      continue;
    }

    path.points = smoothHomotopyPath(
        path.points,
        relevant_obstacle_refs,
        esdf_query,
        d_safe,
        prm_params.m_check,
        prm_params.smooth_spacing,
        prm_params.smooth_passes,
        prm_params.smooth_strength);
    if (path.points.size() < 2 || hasSelfIntersection(path.points)) {
      continue;
    }

    path.length = polylineLength(path.points);
    path.homotopy_signature = signatureKey(
        saturateSignature(computeRawSignatureVector(path.points, relevant_obstacle_refs)));

    auto existing = best_by_signature.find(path.homotopy_signature);
    if (existing == best_by_signature.end() || path.length + 1e-6 < existing->second.length) {
      best_by_signature[path.homotopy_signature] = std::move(path);
    }
  }

  std::vector<GeometricPath> deduplicated_paths;
  deduplicated_paths.reserve(best_by_signature.size());
  for (auto& [signature, path] : best_by_signature) {
    deduplicated_paths.push_back(std::move(path));
  }

  std::sort(
      deduplicated_paths.begin(),
      deduplicated_paths.end(),
      [](const GeometricPath& lhs, const GeometricPath& rhs) {
        return lhs.length < rhs.length;
      });

  if (deduplicated_paths.empty()) {
    last_debug_stats_.unique_homotopy_paths = 0;
    last_debug_stats_.filtered_paths = 0;
    return {};
  }

  const double best_length = deduplicated_paths.front().length;
  const double max_length = std::max(1.0, prm_params.path_length_ratio_limit) * best_length + 1e-6;

  std::vector<GeometricPath> filtered_paths;
  filtered_paths.reserve(deduplicated_paths.size());
  for (auto& path : deduplicated_paths) {
    if (path.length <= max_length) {
      filtered_paths.push_back(std::move(path));
    }
  }

  if (static_cast<int>(filtered_paths.size()) > prm_params.k_max_paths) {
    filtered_paths.resize(static_cast<size_t>(prm_params.k_max_paths));
  }

  last_debug_stats_.unique_homotopy_paths = static_cast<int>(deduplicated_paths.size());
  last_debug_stats_.filtered_paths = static_cast<int>(filtered_paths.size());
  return filtered_paths;
}

}  // namespace mocha
