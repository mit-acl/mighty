#include <mighty/obstacle_tracker_node.hpp>
#include <limits>

// Adaptive EKF Prediction Step for 3D
void ekf_predict(EKFState &ekf_state, double dt)
{
    static Eigen::MatrixXd F = Eigen::MatrixXd::Identity(9, 9);
    F(0, 3) = dt;
    F(1, 4) = dt;
    F(2, 5) = dt;
    F(0, 6) = 0.5 * dt * dt;
    F(1, 7) = 0.5 * dt * dt;
    F(2, 8) = 0.5 * dt * dt;
    F(3, 6) = dt;
    F(4, 7) = dt;
    F(5, 8) = dt;

    ekf_state.x = F * ekf_state.x; // Predict state
    ekf_state.P = F * ekf_state.P.selfadjointView<Eigen::Lower>() * F.transpose() + ekf_state.Q;
}

// Adaptive EKF Update Step for 3D
void aekf_update(EKFState &ekf_state, const Eigen::VectorXd &z, double alpha, double time_updated, const Eigen::Vector3d &bbox, bool use_adaptive_kf)
{
    Eigen::MatrixXd H(3, 9); // Measurement matrix (we only measure position [x, y, z])
    H << 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0;

    Eigen::VectorXd d = z - H * ekf_state.x;                           // Measurement residual
    Eigen::MatrixXd S = H * ekf_state.P * H.transpose() + ekf_state.R; // Residual covariance
    Eigen::MatrixXd K = ekf_state.P * H.transpose() * S.inverse();     // Kalman gain

    // Update state
    ekf_state.x = ekf_state.x + K * d;

    // residual
    Eigen::VectorXd epsilon = z - H * ekf_state.x;

    // Adaptive Kalman Filter
    if (use_adaptive_kf)
    {
        ekf_state.R = alpha * ekf_state.R + (1 - alpha) * (epsilon * epsilon.transpose() + H * ekf_state.P * H.transpose());
        ekf_state.Q = alpha * ekf_state.Q + (1 - alpha) * (K * d * d.transpose() * K.transpose());
    }
    else
    {
        ekf_state.R = Eigen::MatrixXd::Identity(3, 3) * 0.1;
        ekf_state.Q = Eigen::MatrixXd::Identity(9, 9) * 0.01;
    }

    // Update covariance
    ekf_state.P = (Eigen::MatrixXd::Identity(9, 9) - K * H) * ekf_state.P; // Update covariance

    // Update time
    ekf_state.time_updated = time_updated;

    // Update bounding box
    ekf_state.bbox = 0.5 * ekf_state.bbox + (1 - 0.5) * bbox;
}

// Associate cluster with the nearest EKF state using Euclidean distance
int associate_cluster_with_ekf(const Eigen::Vector3d &cluster_centroid, const std::vector<EKFState> &ekf_states, double cluster_tolerance)
{
    double min_distance = cluster_tolerance; // Minimum distance to associate
    int closest_ekf_idx = -1;

    for (int i = 0; i < ekf_states.size(); ++i)
    {
        const Eigen::Vector3d ekf_position(ekf_states[i].x[0], ekf_states[i].x[1], ekf_states[i].x[2]);
        double distance = (cluster_centroid - ekf_position).norm(); // Euclidean distance

        if (distance < min_distance)
        {
            min_distance = distance;
            closest_ekf_idx = i;
        }
    }

    // Return index of the closest EKF state
    return closest_ekf_idx;
}

// ObstacleTrackerNode constructor
ObstacleTrackerNode::ObstacleTrackerNode() : Node("obstacle_tracker_node")
{

    RCLCPP_INFO(this->get_logger(), "Obstacle Tracker Node Started");

    // Declare and set parameters
    declareAndsetParameters();

    // Subscribe to PointCloud2 topic
    sub_pointcloud_ = this->create_subscription<sensor_msgs::msg::PointCloud2>(
        "point_cloud", 10, std::bind(&ObstacleTrackerNode::pointcloudCallback, this, std::placeholders::_1));

    // Publish clusters and predicted positions
    pub_markers_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("tracked_obstacles", 10);
    pub_bboxes_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("cluster_bounding_boxes", 10);

    // Publish the box's location and associated uncertainty as the sphere size
    pub_unc_sphere_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("uncertainty_spheres", 10);

    // Publish predicted trajectory
    pub_predicted_traj_ = this->create_publisher<dynus_interfaces::msg::DynTraj>("predicted_trajs", 10);

    // Initialize the tf2 buffer and listener
    tf2_buffer_ = std::make_shared<tf2_ros::Buffer>(this->get_clock());
    tf2_listener_ = std::make_shared<tf2_ros::TransformListener>(*tf2_buffer_);

    RCLCPP_INFO(this->get_logger(), "Obstacle Tracker Node Initialized");
}

// Declare and set parameters
void ObstacleTrackerNode::declareAndsetParameters()
{

    // Declare parameters
    this->declare_parameter("visual_level", 1);
    this->declare_parameter("use_adaptive_kf", true);
    this->declare_parameter("adaptive_kf_alpha", 0.98);
    this->declare_parameter("adaptive_kf_dt", 0.1);
    this->declare_parameter("cluster_tolerance", 2.0);
    this->declare_parameter("min_cluster_size", 10);
    this->declare_parameter("max_cluster_size", 2000);
    this->declare_parameter("prediction_horizon", 2.0);
    this->declare_parameter("prediction_dt", 0.1);
    this->declare_parameter("time_to_delete_old_obstacles", 10.0);
    this->declare_parameter("cluster_bbox_cutoff_size", 5.0);
    this->declare_parameter("use_life_time_for_box_visualization", false);
    this->declare_parameter("box_visualization_duration", 3.0);
    this->declare_parameter("mighty_map_res", 0.5);
    this->declare_parameter("velocity_threshold", 0.1);
    this->declare_parameter("acceleration_threshold", 0.1);
    this->declare_parameter("use_hardware", false);

    // Set parameters
    visual_level_ = this->get_parameter("visual_level").as_int();
    use_adaptive_kf_ = this->get_parameter("use_adaptive_kf").as_bool();
    adaptive_kf_alpha_ = this->get_parameter("adaptive_kf_alpha").as_double();
    adaptive_kf_dt_ = this->get_parameter("adaptive_kf_dt").as_double();
    cluster_tolerance_ = this->get_parameter("cluster_tolerance").as_double();
    min_cluster_size_ = this->get_parameter("min_cluster_size").as_int();
    max_cluster_size_ = this->get_parameter("max_cluster_size").as_int();
    prediction_horizon_ = this->get_parameter("prediction_horizon").as_double();
    prediction_dt_ = this->get_parameter("prediction_dt").as_double();
    time_to_delete_old_obstacles_ = this->get_parameter("time_to_delete_old_obstacles").as_double();
    cluster_bbox_cutoff_size_ = this->get_parameter("cluster_bbox_cutoff_size").as_double();
    use_life_time_for_box_visualization_ = this->get_parameter("use_life_time_for_box_visualization").as_bool();
    box_visualization_duration_ = this->get_parameter("box_visualization_duration").as_double();
    mighty_map_res_ = this->get_parameter("mighty_map_res").as_double();
    velocity_threshold_ = this->get_parameter("velocity_threshold").as_double();
    acceleration_threshold_ = this->get_parameter("acceleration_threshold").as_double();
    use_hardware_ = this->get_parameter("use_hardware").as_bool();

    // Print the parameters
    // RCLCPP_INFO(this->get_logger(), "visual_level: %d", visual_level_);
    // RCLCPP_INFO(this->get_logger(), "adaptive_kf_alpha: %f", adaptive_kf_alpha_);
    // RCLCPP_INFO(this->get_logger(), "adaptive_kf_dt: %f", adaptive_kf_dt_);
    // RCLCPP_INFO(this->get_logger(), "cluster_tolerance: %f", cluster_tolerance_);
    // RCLCPP_INFO(this->get_logger(), "min_cluster_size: %d", min_cluster_size_);
    // RCLCPP_INFO(this->get_logger(), "max_cluster_size: %d", max_cluster_size_);
    // RCLCPP_INFO(this->get_logger(), "prediction_horizon: %f", prediction_horizon_);
    // RCLCPP_INFO(this->get_logger(), "prediction_dt: %f", prediction_dt_);
    // RCLCPP_INFO(this->get_logger(), "time_to_delete_old_obstacles: %f", time_to_delete_old_obstacles_);
    // RCLCPP_INFO(this->get_logger(), "cluster_bbox_cutoff_size: %f", cluster_bbox_cutoff_size_);
    // RCLCPP_INFO(this->get_logger(), "use_life_time_for_box_visualization: %d", use_life_time_for_box_visualization_);
    // RCLCPP_INFO(this->get_logger(), "box_visualization_duration: %f", box_visualization_duration_);
    // RCLCPP_INFO(this->get_logger(), "mighty_map_res: %f", mighty_map_res_);
    // RCLCPP_INFO(this->get_logger(), "velocity_threshold: %f", velocity_threshold_);
    // RCLCPP_INFO(this->get_logger(), "acceleration_threshold: %f", acceleration_threshold_);
    // RCLCPP_INFO(this->get_logger(), "use_hardware: %d", use_hardware_);
}

// PointCloud2 callback function
void ObstacleTrackerNode::pointcloudCallback(const sensor_msgs::msg::PointCloud2::SharedPtr msg)
{

    // Receive the point cloud message
    // printf("Received point cloud with %zu points\n", msg->width * msg->height);

    // Convert PointCloud2 to PCL format
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::fromROSMsg(*msg, *cloud);

    // Remove NaN values from the cloud
    std::vector<int> indices;
    pcl::removeNaNFromPointCloud(*cloud, *cloud, indices);

    // Voxel grid filtering to downsample the cloud
    pcl::VoxelGrid<pcl::PointXYZ> vg;
    vg.setInputCloud(cloud);
    // vg.setLeafSize(mighty_map_res_, mighty_map_res_, mighty_map_res_);
    vg.setLeafSize(0.2, 0.2, 0.2);
    vg.filter(*cloud);

    // Check if the cloud is empty
    if (cloud->empty())
    {
        RCLCPP_WARN(this->get_logger(), "Received empty point cloud!");
        return;
    }

    // For hardware, we need to remove slash
    std::string targ_frame_id = msg->header.frame_id; 
    if (use_hardware_ && !targ_frame_id.empty() && targ_frame_id[0] == '/')
        targ_frame_id.erase(0, 1);

    // Transform the cloud to the map frame
    geometry_msgs::msg::TransformStamped transform_stamped;
    try
    {
        transform_stamped = tf2_buffer_->lookupTransform(frame_id_,
                                                         targ_frame_id,
                                                         msg->header.stamp,
                                                         rclcpp::Duration::from_seconds(10.0));
    }
    catch (tf2::TransformException &ex)
    {
        RCLCPP_WARN(this->get_logger(), "Transform error: %s", ex.what());
        return;
    }

    // Transform the point cloud to the map frame
    Eigen::Affine3d w_T_b = tf2::transformToEigen(transform_stamped);
    pcl::transformPointCloud(*cloud, *cloud, w_T_b);

    pcl::PassThrough<pcl::PointXYZ> pass_z;
    pass_z.setFilterFieldName("z");
    pass_z.setFilterLimits(0.5, 6.0);
    pass_z.setInputCloud(cloud);
    pass_z.filter(*cloud);

    // Euclidean Cluster Extraction
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
    tree->setInputCloud(cloud);

    std::vector<pcl::PointIndices> cluster_indices;
    pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
    ec.setClusterTolerance(cluster_tolerance_); // Cluster tolerance (distance)
    ec.setMinClusterSize(min_cluster_size_);    // Minimum number of points per cluster
    ec.setMaxClusterSize(max_cluster_size_);    // Maximum number of points per cluster
    ec.setSearchMethod(tree);
    ec.setInputCloud(cloud);
    ec.extract(cluster_indices);

    // Visualize the clustered point cloud
    std::vector<Eigen::Vector3d> cluster_centroids;
    std::vector<Eigen::Vector3d> cluster_bboxes;
    getCentroidsAndSizesOfClusters(cloud, cluster_indices, cluster_centroids, cluster_bboxes);

    std::vector<Cluster> clusters;

    // Delete old EKF states that have not been updated for a long time
    deleteOldEKFstates();

    for (size_t i = 0; i < cluster_indices.size(); ++i)
    {

        auto &indices = cluster_indices[i];

        // Get the cluster's centroid as the measurement
        Eigen::Vector3d centroid = cluster_centroids[i];

        // Get the cluster's bounding box size
        Eigen::Vector3d bbox = cluster_bboxes[i];

        // Check bbox size -> if too large, ignore the cluster (probably a static object)
        if (bbox.norm() > cluster_bbox_cutoff_size_)
        {
            continue;
        }

        // Find the closest EKF state (data association)
        int closest_ekf_idx = associate_cluster_with_ekf(centroid, ekf_states_, cluster_tolerance_);

        // Initialize a new cluster
        Cluster cluster;

        if (closest_ekf_idx >= 0)
        {
            // Update the existing EKF state
            ekf_predict(ekf_states_[closest_ekf_idx], adaptive_kf_dt_);                                                             // EKF Prediction step
            aekf_update(ekf_states_[closest_ekf_idx], centroid, adaptive_kf_alpha_, this->now().seconds(), bbox, use_adaptive_kf_); // EKF Update step
            cluster.setEKFStateAndCentroid(ekf_states_[closest_ekf_idx], centroid);
        }
        else
        {
            // No match found, add a new EKF state
            // Calculate or use averaged Q and R from existing EKF states
            Eigen::MatrixXd Q_avg, R_avg;
            calculateAverageQandR(Q_avg, R_avg);
            EKFState new_state(9, Q_avg, R_avg, this->now().seconds(), bbox, ekf_state_id_++);
            new_state.x.head(3) = centroid; // Initialize state with the centroid
            ekf_states_.push_back(new_state);
            cluster.setEKFStateAndCentroid(new_state, centroid);
        }

        // Add the cluster to the vector
        clusters.push_back(cluster);
    }

    // RCLCPP_INFO(this->get_logger(), "Number of clusters: %d", clusters.size());

    // Publish boxes for visualization
    if (visual_level_ >= 0)
        publishBoxes(clusters);

    // Publish predicted positions (using constant acceleration)
    publishPredictions(clusters);

    // Remove static obstacles based on low velocity
    // filterStaticObstacles();
}

// Function to delete old EKF states that have not been updated for a long time
void ObstacleTrackerNode::deleteOldEKFstates()
{

    // Get the current time
    double current_time = this->now().seconds();

    // Iterate through the EKF states and delete the old ones
    for (size_t i = 0; i < ekf_states_.size(); ++i)
    {
        if (current_time - ekf_states_[i].time_updated > time_to_delete_old_obstacles_)
        {
            ekf_states_.erase(ekf_states_.begin() + i);
        }
    }
}

// Function to calculate the average of Q and R across all EKF states
void ObstacleTrackerNode::calculateAverageQandR(Eigen::MatrixXd &Q_avg, Eigen::MatrixXd &R_avg)
{
    Q_avg = Eigen::MatrixXd::Zero(9, 9);
    R_avg = Eigen::MatrixXd::Zero(3, 3);

    // TODO: we can implemente weighted average based on the time since the last update (old obstacles have less weight)
    if (!ekf_states_.empty())
    {
        for (const auto &ekf_state : ekf_states_)
        {
            Q_avg += ekf_state.Q;
            R_avg += ekf_state.R;
        }
        Q_avg /= ekf_states_.size();
        R_avg /= ekf_states_.size();
    }
    else
    {
        // If there are no EKF states, initialize Q and R to default values
        Q_avg = Eigen::MatrixXd::Identity(9, 9) * 0.01;
        R_avg = Eigen::MatrixXd::Identity(3, 3) * 0.1;
    }
}

void ObstacleTrackerNode::getCentroidsAndSizesOfClusters(const pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, const std::vector<pcl::PointIndices> &cluster_indices, std::vector<Eigen::Vector3d> &cluster_centroids, std::vector<Eigen::Vector3d> &cluster_bboxes)
{

    // Create and add new markers & Get the min/max values for each cluster
    for (size_t i = 0; i < cluster_indices.size(); ++i)
    {
        auto &indices = cluster_indices[i];

        // Use Eigen vectorization
        Eigen::Vector3d min = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
        Eigen::Vector3d max = Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest());

        for (const auto &idx : indices.indices)
        {
            const auto &point = cloud->points[idx];
            min = min.cwiseMin(Eigen::Vector3d(point.x, point.y, point.z));
            max = max.cwiseMax(Eigen::Vector3d(point.x, point.y, point.z));
        }

        Eigen::Vector3d centroid = (min + max) / 2.0;
        Eigen::Vector3d bbox = max - min;

        // Add centroid and size to the vectors
        cluster_centroids.emplace_back(centroid);
        cluster_bboxes.emplace_back(bbox);
    }
}

void ObstacleTrackerNode::publishBoxes(const std::vector<Cluster> &clusters)
{

    // printf("publishBoxes called with %zu clusters\n", clusters.size());

    visualization_msgs::msg::MarkerArray cluster_markers;
    visualization_msgs::msg::MarkerArray unc_sphere_markers;

    // Set the max scale to avoid very large spheres
    double max_scale = 2.5;

    // Create and add new markers & Get the min/max values for each cluster
    for (auto &cluster : clusters)
    {
        // Create a CUBE marker to visualize the bounding box
        visualization_msgs::msg::Marker marker;
        if (use_life_time_for_box_visualization_)
            marker.lifetime = rclcpp::Duration::from_seconds(box_visualization_duration_);
        marker.header.frame_id = frame_id_;
        marker.header.stamp = this->now(); // Ensure timestamp consistency
        marker.ns = "cluster_bounding_box";
        marker.id = marker_id_++; // Ensure unique IDs for new markers
        marker.type = visualization_msgs::msg::Marker::CUBE;
        marker.action = visualization_msgs::msg::Marker::ADD;

        // Set the position to the center of the bounding box
        marker.pose.position.x = cluster.centroid[0];
        marker.pose.position.y = cluster.centroid[1];
        marker.pose.position.z = cluster.centroid[2];

        // Set the scale to the size of the bounding box
        marker.scale.x = cluster.ekf_state.bbox[0];
        marker.scale.y = cluster.ekf_state.bbox[1];
        marker.scale.z = cluster.ekf_state.bbox[2];

        // Set the color
        marker.color = cluster.ekf_state.color;

        // Set alpha (transparency)
        marker.color.a = 1.0;

        // Add the bounding box marker to the marker array
        cluster_markers.markers.push_back(marker);

        // Create a SPHERE marker to visualize the uncertainty
        visualization_msgs::msg::Marker unc_sphere_marker;
        // if (use_life_time_for_box_visualization_) unc_sphere_marker.lifetime = rclcpp::Duration::from_seconds(box_visualization_duration_);
        unc_sphere_marker.header.frame_id = frame_id_;
        unc_sphere_marker.header.stamp = this->now(); // Ensure timestamp consistency
        unc_sphere_marker.ns = "uncertainty_sphere";
        unc_sphere_marker.id = 0; // Not unique IDs so that only one sphere is visualized
        unc_sphere_marker.type = visualization_msgs::msg::Marker::SPHERE;
        unc_sphere_marker.action = visualization_msgs::msg::Marker::ADD;

        // Set the position to the center of the bounding box
        unc_sphere_marker.pose.position.x = cluster.centroid[0];
        unc_sphere_marker.pose.position.y = cluster.centroid[1];
        unc_sphere_marker.pose.position.z = cluster.centroid[2];

        // Set the scale to the size of uncertainty
        unc_sphere_marker.scale.x = std::min(cluster.ekf_state.P.diagonal()[0] * 2e3, max_scale);
        unc_sphere_marker.scale.y = std::min(cluster.ekf_state.P.diagonal()[1] * 2e3, max_scale);
        unc_sphere_marker.scale.z = std::min(cluster.ekf_state.P.diagonal()[2] * 2e3, max_scale);

        // Set the color
        unc_sphere_marker.color = cluster.ekf_state.color;

        // If you want to set a specific color
        // unc_sphere_marker.color.r = 255.0 / 255.0;
        // unc_sphere_marker.color.g = 51.0 / 255.0;
        // unc_sphere_marker.color.b = 153.0 / 255.0;

        // Set alpha (transparency)
        unc_sphere_marker.color.a = 0.6;

        // Add the uncertainty sphere marker to the marker array
        unc_sphere_markers.markers.push_back(unc_sphere_marker);
    }

    // Step 3: Publish the marker array to visualize clusters
    pub_bboxes_->publish(cluster_markers);

    // Step 4: Publish the marker array to visualize uncertainty spheres
    pub_unc_sphere_->publish(unc_sphere_markers);
}

// Visualize the predictions over a time horizon, updating both position and velocity
void ObstacleTrackerNode::publishPredictions(const std::vector<Cluster> &clusters)
{

    visualization_msgs::msg::MarkerArray markers;
    int id = 0;
    int num_steps = static_cast<int>(prediction_horizon_ / prediction_dt_); // Number of steps

    // Initialize knots and positions
    std::vector<double> t_values;
    std::vector<double> x_values;
    std::vector<double> y_values;
    std::vector<double> z_values;

    for (int i = 0; i < clusters.size(); ++i)
    {

        // Initialize position, velocity, and acceleration from EKF state
        Eigen::Vector3d current_position(clusters[i].ekf_state.x[0], clusters[i].ekf_state.x[1], clusters[i].ekf_state.x[2]);
        Eigen::Vector3d current_velocity(clusters[i].ekf_state.x[3], clusters[i].ekf_state.x[4], clusters[i].ekf_state.x[5]);
        Eigen::Vector3d acceleration(clusters[i].ekf_state.x[6], clusters[i].ekf_state.x[7], clusters[i].ekf_state.x[8]);

        // Avoid high acceleration values
        for (int j = 0; j < 3; ++j)
        {
            if (acceleration[j] > 2.0)
                acceleration[j] = 2.0;
            if (acceleration[j] < -2.0)
                acceleration[j] = -2.0;
        }

        // Store the initial position and time
        t_values.push_back(0.0);
        x_values.push_back(current_position[0]);
        y_values.push_back(current_position[1]);
        z_values.push_back(current_position[2]);

        for (int step = 0; step < num_steps; ++step)
        {

            for (int j = 0; j < 3; ++j)
            {
                // Avoid high velocity values
                if (current_velocity[j] > 0.5)
                    current_velocity[j] = 0.5;
                if (current_velocity[j] < -0.5)
                    current_velocity[j] = -0.5;
            
                // Avoid high acceleration values
                if (acceleration[j] > 0.8)
                    acceleration[j] = 0.8;
                if (acceleration[j] < -0.8)
                    acceleration[j] = -0.8;
            }

            // Calculate time for this step
            double t = step * prediction_dt_;

            // Update velocity: v(t) = v_0 + a * t
            Eigen::Vector3d future_velocity = current_velocity + acceleration * t;

            // Predict future position: p(t) = p_0 + v_0 * t + 0.5 * a * t^2
            Eigen::Vector3d future_position;
            future_position[0] = current_position[0] + current_velocity[0] * t + 0.5 * acceleration[0] * t * t;
            future_position[1] = current_position[1] + current_velocity[1] * t + 0.5 * acceleration[1] * t * t;
            future_position[2] = current_position[2] + current_velocity[2] * t + 0.5 * acceleration[2] * t * t;

            // Store the time and position values
            t_values.push_back(t);
            x_values.push_back(future_position[0]);
            y_values.push_back(future_position[1]);
            z_values.push_back(future_position[2]);

            // Create a marker to visualize the predicted position at this time step
            visualization_msgs::msg::Marker marker;
            marker.header.frame_id = frame_id_; // Ensure it's the correct frame
            marker.id = id++;
            marker.type = visualization_msgs::msg::Marker::ARROW;
            marker.action = visualization_msgs::msg::Marker::ADD;

            // Set arrow start (current position) and end (future position)
            geometry_msgs::msg::Point start, end;
            start.x = current_position[0];
            start.y = current_position[1];
            start.z = current_position[2];
            end.x = future_position[0];
            end.y = future_position[1];
            end.z = future_position[2];

            // Set arrow start and end points
            marker.points.push_back(start);
            marker.points.push_back(end);

            // Set scale (arrow width and length)
            marker.scale.x = 0.1;
            marker.scale.y = 0.2;

            // Set color (you can gradually fade it based on time step)
            marker.color.r = clusters[i].ekf_state.color.r;
            marker.color.g = clusters[i].ekf_state.color.g;
            marker.color.b = clusters[i].ekf_state.color.b;
            marker.color.a = clusters[i].ekf_state.color.a;

            // Add this marker to the marker array
            markers.markers.push_back(marker);

            // Update the current position and velocity for the next step
            current_position = future_position;
            current_velocity = future_velocity;
        }

        // Check if x_values, y_values, z_values changed much (especially in the beginning, the predicted trajectories can be very short and hard to fit)
        // If they are not changing much, we can skip the polynomial fitting
        double x_diff = abs(x_values.front() - x_values.back());
        double y_diff = abs(y_values.front() - y_values.back());
        double z_diff = abs(z_values.front() - z_values.back());

        double cutoff_length_threshold = 0.1; // TODO: make this a parameter?

        // If the predicted trajectory is too short, skip the polynomial fitting
        if (x_diff < cutoff_length_threshold && y_diff < cutoff_length_threshold && z_diff < cutoff_length_threshold)
        {
            t_values.clear();
            x_values.clear();
            y_values.clear();
            z_values.clear();
            continue;
        }

        // Fit a polynomial to the predicted positions
        Eigen::VectorXd beta_x = polyfit(t_values, x_values, degree_for_pwp_);
        Eigen::VectorXd beta_y = polyfit(t_values, y_values, degree_for_pwp_);
        Eigen::VectorXd beta_z = polyfit(t_values, z_values, degree_for_pwp_);

        // Calculate variance of the residuals
        double variance_x = calculateVariance(t_values, x_values, beta_x, degree_for_pwp_);
        double variance_y = calculateVariance(t_values, y_values, beta_y, degree_for_pwp_);
        double variance_z = calculateVariance(t_values, z_values, beta_z, degree_for_pwp_);

        // Convert t_values, beta_x, beta_y, beta_z to PieceWisePol
        PieceWisePol pwp;
        double current_time = this->now().seconds();
        pwp.times.push_back(current_time);
        pwp.times.push_back(current_time + prediction_horizon_);
        pwp.coeff_x.push_back({beta_x(0), beta_x(1), beta_x(2), beta_x(3)});
        pwp.coeff_y.push_back({beta_y(0), beta_y(1), beta_y(2), beta_y(3)});
        pwp.coeff_z.push_back({beta_z(0), beta_z(1), beta_z(2), beta_z(3)});

        // Fit a quintic polynomial to the predicted positions
        Eigen::VectorXd beta_x_quintic = polyfit(t_values, x_values, degree_for_poly_);
        Eigen::VectorXd beta_y_quintic = polyfit(t_values, y_values, degree_for_poly_);
        Eigen::VectorXd beta_z_quintic = polyfit(t_values, z_values, degree_for_poly_);

        // Publish DynTraj message with the predicted trajectory
        dynus_interfaces::msg::DynTraj msg;
        msg.header.stamp = this->now();
        msg.header.frame_id = frame_id_;
        msg.id = clusters[i].ekf_state.id;
        msg.bbox.push_back(clusters[i].ekf_state.bbox.x());
        msg.bbox.push_back(clusters[i].ekf_state.bbox.y());
        msg.bbox.push_back(clusters[i].ekf_state.bbox.z());
        msg.pwp = mighty_utils::convertPwp2PwpMsg(pwp);
        msg.ekf_cov_p.push_back(clusters[i].ekf_state.P(0, 0));
        msg.ekf_cov_p.push_back(clusters[i].ekf_state.P(1, 1));
        msg.ekf_cov_p.push_back(clusters[i].ekf_state.P(2, 2));
        msg.ekf_cov_q.push_back(clusters[i].ekf_state.Q(0, 0));
        msg.ekf_cov_q.push_back(clusters[i].ekf_state.Q(1, 1));
        msg.ekf_cov_q.push_back(clusters[i].ekf_state.Q(2, 2));
        msg.ekf_cov_r.push_back(clusters[i].ekf_state.R(0, 0));
        msg.ekf_cov_r.push_back(clusters[i].ekf_state.R(1, 1));
        msg.ekf_cov_r.push_back(clusters[i].ekf_state.R(2, 2));
        msg.poly_cov.push_back(variance_x);
        msg.poly_cov.push_back(variance_y);
        msg.poly_cov.push_back(variance_z);

        // coefficients for quintic polynomial
        msg.poly_coeffs_x.clear();
        msg.poly_coeffs_y.clear();
        msg.poly_coeffs_z.clear();
        for (int j = 0; j < degree_for_poly_ + 1; ++j)
        {
            msg.poly_coeffs_x.push_back(beta_x_quintic(j));
            msg.poly_coeffs_y.push_back(beta_y_quintic(j));
            msg.poly_coeffs_z.push_back(beta_z_quintic(j));
        }

        // Set the start and end times for the trajectory
        msg.poly_start_time = current_time;
        msg.poly_end_time = current_time + prediction_horizon_;

        msg.is_agent = false;
        pub_predicted_traj_->publish(msg);

        // Clear the vectors for the next EKF state
        t_values.clear();
        x_values.clear();
        y_values.clear();
        z_values.clear();
    }

    // Publish the marker array with predicted trajectories
    if (visual_level_ >= 0)
        pub_markers_->publish(markers);
}

Eigen::VectorXd ObstacleTrackerNode::polyfit(const std::vector<double> &t, const std::vector<double> &y, int degree)
{

    // Number of data points
    int n = t.size();

    // Construct the Vandermonde matrix for polynomial fitting
    Eigen::MatrixXd X(n, degree + 1);
    Eigen::VectorXd Y(n);

    for (int i = 0; i < n; ++i)
    {
        Y(i) = y[i];
        for (int j = 0; j <= degree; ++j)
        {
            X(i, j) = std::pow(t[i], degree - j); // Note: the coefficients are stored as [a b c d] for a*t^3 + b*t^2 + c*t + d so we need to flip it instead of std::pow(t[i], j)
        }
    }

    // Solve the normal equations: X^T * X * beta = X^T * Y
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * Y);

    return beta;
}

double ObstacleTrackerNode::calculateVariance(const std::vector<double> &t, const std::vector<double> &y, const Eigen::VectorXd &beta, int degree)
{
    // Calculate residuals and estimate variance
    int n = t.size();
    double residual_sum = 0.0;

    for (int i = 0; i < n; ++i)
    {
        double fitted_value = 0.0;
        for (int j = 0; j <= degree; ++j)
        {
            fitted_value += beta(j) * std::pow(t[i], j);
        }
        double residual = y[i] - fitted_value;
        residual_sum += residual * residual;
    }

    return residual_sum / (n - degree - 1); // variance = sum(residuals^2) / (n - p - 1)
}

// Function to filter out static obstacles based on low velocity
void ObstacleTrackerNode::filterStaticObstacles()
{
    for (auto it = ekf_states_.begin(); it != ekf_states_.end();)
    {
        Eigen::Vector3d velocity = it->x.segment<3>(3); // Extract velocity from the state vector (elements 3,4,5)
        Eigen::Vector3d acceleration = it->x.segment<3>(6); // Extract acceleration from the state vector (elements 6,7,8)
        if (velocity.norm() < velocity_threshold_ && acceleration.norm() < acceleration_threshold_)
        {
            // If the velocity is below the threshold, consider the obstacle static and remove it
            it = ekf_states_.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

int main(int argc, char **argv)
{

    rclcpp::init(argc, argv);

    // Initialize multi-threaded executor
    rclcpp::executors::MultiThreadedExecutor executor;

    // Create an instance of ObstacleTrackerNode
    auto obstacle_tracker_node = std::make_shared<ObstacleTrackerNode>();
    executor.add_node(obstacle_tracker_node);

    // Spin the executor
    executor.spin();

    rclcpp::shutdown();
    return 0;
}
