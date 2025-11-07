/**
 * @file graph_search.h
 * @brief backend of graph search, implementation of A* and JPS
 */

#ifndef DGP_GRAPH_SEARCH_H
#define DGP_GRAPH_SEARCH_H

#include <boost/heap/d_ary_heap.hpp> // boost::heap::d_ary_heap
#include <memory>                    // std::shared_ptr
#include <limits>                    // std::numeric_limits
#include <vector>                    // std::vector
#include <unordered_map>             // std::unordered_map
#include <mutex>                     // std::mutex
#include <dgp/map_util.hpp>          // mighty::MapUtil
#include <dgp/data_type.hpp>
#include <timer.hpp>
#include <dgp/data_type.hpp>
#include <mighty/utils.hpp>
#include <mighty/mighty_type.hpp>

// TODO: ROS dependency (not ideal)
#include <rclcpp/rclcpp.hpp>

namespace mighty
{
  /// Heap element comparison
  template <class T>
  struct compare_state
  {
    bool operator()(T a1, T a2) const
    {
      double f1 = a1->g + a1->h;
      double f2 = a2->g + a2->h;
      if ((f1 >= f2 - 0.000001) && (f1 <= f2 + 0.000001))
        return a1->g < a2->g; // if equal compare gvals
      return f1 > f2;
    }
  };

  /// Define priority queue
  struct State; // forward declaration
  /// State pointer
  using StatePtr = std::shared_ptr<State>;
  using priorityQueue = boost::heap::d_ary_heap<StatePtr, boost::heap::mutable_<true>,
                                                boost::heap::arity<2>, boost::heap::compare<compare_state<StatePtr>>>;

  /// Node of the graph in graph search
  struct State
  {
    /// ID
    int id;
    /// Coord
    int x, y, z = 0;
    /// direction
    int dx, dy, dz; // discrete coordinates of this node
    /// id of predicessors
    int parentId = -1;

    /// pointer to heap location
    priorityQueue::handle_type heapkey;

    /// g cost
    double g = std::numeric_limits<double>::infinity();
    /// heuristic cost
    double h;
    /// if has been opened
    bool opened = false;
    /// if has been closed
    bool closed = false;

    /// travel_time used for dynamic A*
    double t = 0.0;

    /// velocity used for dynamic A*
    Eigen::Vector3d vel = Eigen::Vector3d::Zero();

    /// 2D constructor
    State(int id, int x, int y, int dx, int dy)
        : id(id), x(x), y(y), dx(dx), dy(dy)
    {
    }

    /// 3D constructor
    State(int id, int x, int y, int z, int dx, int dy, int dz)
        : id(id), x(x), y(y), z(z), dx(dx), dy(dy), dz(dz)
    {
    }

    /// 3D constructor with g value
    State(int id, int x, int y, int z, int dx, int dy, int dz, double g)
        : id(id), x(x), y(y), z(z), dx(dx), dy(dy), dz(dz), g(g)
    {
    }
  };

  /// Search and prune neighbors for JPS 2D
  struct JPS2DNeib
  {
    // for each (dx,dy) these contain:
    //    ns: neighbors that are always added
    //    f1: forced neighbors to check
    //    f2: neighbors to add if f1 is forced
    int ns[9][2][8];
    int f1[9][2][2];
    int f2[9][2][2];
    // nsz contains the number of neighbors for the four different types of moves:
    // no move (norm 0):        8 neighbors always added
    //                          0 forced neighbors to check (never happens)
    //                          0 neighbors to add if forced (never happens)
    // straight (norm 1):       1 neighbor always added
    //                          2 forced neighbors to check
    //                          2 neighbors to add if forced
    // diagonal (norm sqrt(2)): 3 neighbors always added
    //                          2 forced neighbors to check
    //                          2 neighbors to add if forced
    static constexpr int nsz[3][2] = {{8, 0}, {1, 2}, {3, 2}};

    void print();
    JPS2DNeib();

  private:
    void Neib(int dx, int dy, int norm1, int dev, int &tx, int &ty);
    void FNeib(int dx, int dy, int norm1, int dev,
               int &fx, int &fy, int &nx, int &ny);
  };

  /// Search and prune neighbors for JPS 3D
  struct JPS3DNeib
  {
    // for each (dx,dy,dz) these contain:
    //    ns: neighbors that are always added
    //    f1: forced neighbors to check
    //    f2: neighbors to add if f1 is forced
    int ns[27][3][26];
    int f1[27][3][12];
    int f2[27][3][12];
    // nsz contains the number of neighbors for the four different types of moves:
    // no move (norm 0):        26 neighbors always added
    //                          0 forced neighbors to check (never happens)
    //                          0 neighbors to add if forced (never happens)
    // straight (norm 1):       1 neighbor always added
    //                          8 forced neighbors to check
    //                          8 neighbors to add if forced
    // diagonal (norm sqrt(2)): 3 neighbors always added
    //                          8 forced neighbors to check
    //                          12 neighbors to add if forced
    // diagonal (norm sqrt(3)): 7 neighbors always added
    //                          6 forced neighbors to check
    //                          12 neighbors to add if forced
    static constexpr int nsz[4][2] = {{26, 0}, {1, 8}, {3, 12}, {7, 12}};
    JPS3DNeib();

  private:
    void Neib(int dx, int dy, int dz, int norm1, int dev, int &tx, int &ty, int &tz);
    void FNeib(int dx, int dy, int dz, int norm1, int dev,
               int &fx, int &fy, int &fz,
               int &nx, int &ny, int &nz);
  };

  /**
   * @brief GraphSearch class
   *
   * Implement A* and Jump Point Search
   */
  class GraphSearch
  {
  public:
    /**
     * @brief 3D graph search constructor
     *
     * @param cMap 1D array stores the occupancy, with the order equal to \f$x + xDim * y + xDim * yDim * z\f$
     * @param map_util map util
     * @param xDim map length
     * @param yDim map width
     * @param zDim map height
     * @param eps weight of heuristic, optional, default as 1
     * @param verbose flag for printing debug info, optional, default as False
     * @param global_planner initial guess planner, optional, default as ""
     * @param res map resolution, optional, default as 0.5
     */
    GraphSearch(const int *cMap, const std::shared_ptr<mighty::VoxelMapUtil> &map_util, int xDim, int yDim, int zDim, double eps, bool verbose, std::string global_planner, double w_unknown, double w_align = 60.0, double decay_len_cells = 20.0, double w_side = 0.2);

    /**
     * @brief start 3D planning thread
     *
     * @param xStart start x coordinate
     * @param yStart start y coordinate
     * @param zStart start z coordinate
     * @param xGoal goal x coordinate
     * @param yGoal goal y coordinate
     * @param zGoal goal z coordinate
     * @param initial_g initial g value
     * @param max_expand maximum number of expansion allowed, optional, default is -1, means no limitation
     */
    bool plan(int xStart, int yStart, int zStart, int xGoal, int yGoal, int zGoal, double initial_g, double &global_planning_time, double &dgp_static_jps_time, double &dgp_check_path_time, double &dgp_dynamic_astar_time, double &dgp_recover_path_time, double current_time, const Vec3f &start_vel, int max_expand = -1, int timeout_duration_ms = 1000);

    /// Get the optimal path
    std::vector<StatePtr> getPath() const;

    /// Get the states in open set
    std::vector<StatePtr> getOpenSet() const;

    /// Get the states in close set
    std::vector<StatePtr> getCloseSet() const;

    /// Get the states in hash map
    std::vector<StatePtr> getAllSet() const;

    /// Set start and goal
    void setStartAndGoal(const Vecf<3> &start, const Vecf<3> &goal);

    /// Set max dynamic constraints
    void setBounds(double max_values[3]);

  private:
    /// Select planner
    bool select_planner(StatePtr &currNode_ptr, int max_expand, int start_id, int goal_id, std::chrono::milliseconds timeout_duration);
    /// Main planning loop for Static JPS
    bool static_jps_plan(StatePtr &currNode_ptr, int max_expand, int start_id, int goal_id, std::chrono::milliseconds timeout_duration);
    /// Main planning loop for Dynamic A*
    bool dynamic_astar_plan(StatePtr &currNode_ptr, int max_expand, int start_id, int goal_id, std::chrono::milliseconds timeout_duration);
    /// Main planning loop for Dynamic JPS++
    bool dgp_plan(StatePtr &currNode_ptr, int max_expand, int start_id, int goal_id, std::chrono::milliseconds timeout_duration);
    /// Get successor function for A*
    void getSucc(const StatePtr &curr, std::vector<int> &succ_ids, std::vector<double> &succ_costs);
    /// Get successor function for JPS
    void getJpsSucc(const StatePtr &curr, std::vector<int> &succ_ids, std::vector<double> &succ_costs);
    /// Recover the optimal path
    std::vector<StatePtr> recoverPath(StatePtr node, int id);

    /// Get subscript
    int coordToId(int x, int y, int z) const;

    /// Check if (x, y, z) is free
    bool isFree(int x, int y, int z) const;

    /// Check if (x, y, z) is occupied
    bool isOccupied(int x, int y, int z) const;

    /// Check if (x, y, z) is unknown
    bool isUnknown(int x, int y, int z) const;

    /// Clculate heuristic
    double getHeur(int x, int y, int z) const;

    /// Determine if (x, y, z) has forced neighbor with direction (dx, dy, dz)
    bool hasForced(int x, int y, int z, int dx, int dy, int dz);

    /// 3D jump, return true iff finding the goal or a jump point
    bool jump(int x, int y, int z, int dx, int dy, int dz, int &new_x, int &new_y, int &new_z);

    /// Compute node time
    double computeNodeTime(const Eigen::Vector3d &pf, const Eigen::Vector3d &p0, const Eigen::Vector3d &v0, StatePtr &currNode_ptr);

    /// Simplify the path
    std::vector<StatePtr> removeLinePts(const std::vector<StatePtr> &path);
    std::vector<StatePtr> removeCornerPts(const std::vector<StatePtr> &path);

    /// Update G values
    void updateGValues();

    inline StatePtr &atHM(int id)
    {
      static StatePtr nullPtr = nullptr;
      if (id < 0 || id >= (int)hm_.size())
        return nullPtr;
      return hm_[id];
    }

    /// cMap pointer
    const int *cMap_;

    int xDim_, yDim_, zDim_;
    double eps_;
    bool verbose_;

    const int val_free_ = 0;
    const int val_occupied_ = 100;
    const int val_unknown_ = -1;

    const double w_unknown_ = 1.5; // cost for unknown cells, in multiples of base cost (1 for orthogonal step)
    // weights are in "cell-length" units to match succ_costs
    double w_align_{60.0};     // strength of alignment penalty (cells)
    double decay_len_cells_{20.0}; // e-folding distance from the start (cells)
    double w_side_{0.2};      // side (handedness) tie-break strength (cells)

    int xGoal_, yGoal_, zGoal_;
    bool use_jps_ = false;

    priorityQueue pq_;
    std::vector<StatePtr> hm_;
    std::vector<bool> seen_;

    std::vector<StatePtr> path_;

    std::vector<std::vector<int>> ns_;
    std::shared_ptr<JPS2DNeib> jn2d_;
    std::shared_ptr<JPS3DNeib> jn3d_;

    // Initial guess planner
    std::string global_planner_;

    // Map util
    std::shared_ptr<mighty::VoxelMapUtil> map_util_;

    // Set start and goal
    Vecf<3> start_;
    Vecf<3> goal_;

    // Dynamic constraints
    double v_max_;
    Eigen::Vector3d v_max_3d_;
    double a_max_;
    Eigen::Vector3d a_max_3d_;
    double j_max_;
    Eigen::Vector3d j_max_3d_;

    // For benchmarking time recording
    double global_planning_time_ = 0.0;
    double dgp_static_jps_time_ = 0.0;
    double dgp_check_path_time_ = 0.0;
    double dgp_dynamic_astar_time_ = 0.0;
    double dgp_recover_path_time_ = 0.0;

    // For dynamic map
    double current_time_ = 0.0;

    // variables
    Vecf<3> start_vel_;

  };
}
#endif