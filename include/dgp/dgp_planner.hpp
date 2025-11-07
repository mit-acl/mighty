/**
 * @file dgp_planner.h
 * @brief DGPPlanner
 */
#ifndef DGP_PLANNER_BASE_H
#define DGP_PLANNER_BASE_H

#include <dgp/data_type.hpp>
#include <dgp/graph_search.hpp>
#include <dgp/map_util.hpp>

/**
 * @brief Abstract base for planning
 */
class DGPPlanner
{
public:
  /**
   * @brief Simple constructor
   * @param verbose enable debug mode
   */
  DGPPlanner(std::string global_planner, bool verbose, double v_max, double a_max, double j_max, int dgp_timeout_duration_ms, double w_unknown, double w_align, double decay_len_cells, double w_side, int los_cells = 3, double min_len = 0.5, double min_turn = 10.0);

  // Set map util for collistion checking
  void setMapUtil(const std::shared_ptr<mighty::MapUtil<3>> &map_util);
  /**
   * @brief Status of the planner
   *
   * 0 --- exit normally;
   * -1 --- no path found;
   * 1, 2 --- start or goal is not free.
   */
  int status();
  // Get the modified path
  vec_Vecf<3> getPath();
  // Get the raw path
  vec_Vecf<3> getRawPath();
  // Set timeout duration
  void setTimeOutDurationMs(int timeout_duration_ms);
  // remove redundant points on the same line
  vec_Vecf<3> removeLinePts(const vec_Vecf<3> &path);
  // Remove some corner waypoints
  vec_Vecf<3> removeCornerPts(const vec_Vecf<3> &path);
  // Clean up path
  void cleanUpPath(vec_Vecf<3> &path);
  // Planning function
  bool plan(const Vecf<3> &start, const Vecf<3> &start_vel, const Vecf<3> &goal, double &final_g, double current_time, decimal_t eps = 1);
  // Get the nodes in open set
  vec_Vecf<3> getOpenSet() const;
  // Get the nodes in close set
  vec_Vecf<3> getCloseSet() const;
  // Get all the nodes
  vec_Vecf<3> getAllSet() const;
  // Get computation time
  double getInitialGuessPlanningTime();
  double getStaticJPSPlanningTime();
  double getCheckPathTime();
  double getDynamicAstarTime();
  double getRecoverPathTime();
  void updateVmax(double v_max);
  vec_Vecf<3> shortCutByLoS(const vec_Vecf<3> &in, int inflate_radius_cells) const;
  bool lineOfSightCapsule(const Vecf<3> &a,
                          const Vecf<3> &b,
                          int inflate_radius_cells) const;
  vec_Vecf<3> collapseShortEdges(const vec_Vecf<3> &in, double min_len) const;
  vec_Vecf<3> angleSpacingFilter(const vec_Vecf<3> &in,
                                 double min_turn_deg,
                                 double min_seg_len) const;

protected:
  // Assume using 3D voxel map for all 2d and 3d planning
  std::shared_ptr<mighty::VoxelMapUtil> map_util_;
  // The planner
  std::shared_ptr<mighty::GraphSearch> graph_search_;
  // Raw path from planner
  vec_Vecf<3> raw_path_;
  // Modified path for future usage
  vec_Vecf<3> path_;
  // Flag indicating the success of planning
  int status_ = 0;
  // Enabled for printing info
  bool planner_verbose_;

  // Parameters
  std::string global_planner_;

  // For benchmarking time recording
  double global_planning_time_ = 0.0;
  double dgp_static_jps_time_ = 0.0;
  double dgp_check_path_time_ = 0.0;
  double dgp_dynamic_astar_time_ = 0.0;
  double dgp_recover_path_time_ = 0.0;

  // time out duration
  int dgp_timeout_duration_ms_ = 1000;

  // max values
  double v_max_;
  double a_max_;
  double j_max_;

  // alignment
  double w_align_ = 60.0;        // weight for alignment with start velocity
  double decay_len_cells_ = 20.0;    // decay length in cells for alignment
  double w_side_ = 0.2;         // weight for side tie-break

  // unknown cell cost weight
  double w_unknown_ = 1.0;

  // LOS post processing
  int los_cells_ = 3;           // number of cells for inflation in LoS
  double min_len_ = 0.5;        // minimum length of edges
  double min_turn_ = 10.0;      // minimum turn angle in degrees

};

#endif
