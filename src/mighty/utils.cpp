/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "mighty/utils.hpp"

namespace mighty_utils
{

  dynus_interfaces::msg::PWPTraj convertPwp2PwpMsg(const PieceWisePol &pwp)
  {
    dynus_interfaces::msg::PWPTraj pwp_msg;

    for (int i = 0; i < pwp.times.size(); i++)
    {
      pwp_msg.times.push_back(pwp.times[i]);
    }

    // push x
    for (auto coeff_x_i : pwp.coeff_x)
    {
      dynus_interfaces::msg::CoeffPoly3 coeff_poly3;
      coeff_poly3.a = coeff_x_i(0);
      coeff_poly3.b = coeff_x_i(1);
      coeff_poly3.c = coeff_x_i(2);
      coeff_poly3.d = coeff_x_i(3);
      pwp_msg.coeff_x.push_back(coeff_poly3);
    }

    // push y
    for (auto coeff_y_i : pwp.coeff_y)
    {
      dynus_interfaces::msg::CoeffPoly3 coeff_poly3;
      coeff_poly3.a = coeff_y_i(0);
      coeff_poly3.b = coeff_y_i(1);
      coeff_poly3.c = coeff_y_i(2);
      coeff_poly3.d = coeff_y_i(3);
      pwp_msg.coeff_y.push_back(coeff_poly3);
    }

    // push z
    for (auto coeff_z_i : pwp.coeff_z)
    {
      dynus_interfaces::msg::CoeffPoly3 coeff_poly3;
      coeff_poly3.a = coeff_z_i(0);
      coeff_poly3.b = coeff_z_i(1);
      coeff_poly3.c = coeff_z_i(2);
      coeff_poly3.d = coeff_z_i(3);
      pwp_msg.coeff_z.push_back(coeff_poly3);
    }

    return pwp_msg;
  }

  PieceWisePol convertPwpMsg2Pwp(const dynus_interfaces::msg::PWPTraj &pwp_msg)
  {
    PieceWisePol pwp;

    if (pwp_msg.coeff_x.size() != pwp_msg.coeff_y.size() || pwp_msg.coeff_x.size() != pwp_msg.coeff_z.size())
    {
      std::cout << " coeff_x,coeff_y,coeff_z of pwp_msg should have the same elements" << std::endl;
      std::cout << " ================================" << std::endl;
      abort();
    }

    for (int i = 0; i < pwp_msg.times.size(); i++)
    {
      pwp.times.push_back(pwp_msg.times[i]);
    }

    for (int i = 0; i < pwp_msg.coeff_x.size(); i++)
    {
      Eigen::Matrix<double, 4, 1> tmp_x, tmp_y, tmp_z;

      tmp_x << pwp_msg.coeff_x[i].a, pwp_msg.coeff_x[i].b, pwp_msg.coeff_x[i].c, pwp_msg.coeff_x[i].d;
      pwp.coeff_x.push_back(tmp_x);

      tmp_y << pwp_msg.coeff_y[i].a, pwp_msg.coeff_y[i].b, pwp_msg.coeff_y[i].c, pwp_msg.coeff_y[i].d;
      pwp.coeff_y.push_back(tmp_y);

      tmp_z << pwp_msg.coeff_z[i].a, pwp_msg.coeff_z[i].b, pwp_msg.coeff_z[i].c, pwp_msg.coeff_z[i].d;
      pwp.coeff_z.push_back(tmp_z);
    }

    return pwp;
  }

  visualization_msgs::msg::MarkerArray convertPwp2ColoredMarkerArray(PieceWisePol &pwp, int samples)
  {

    // Initialize the marker array
    visualization_msgs::msg::MarkerArray marker_array;

    // Set the initial and final time and get the delta time
    double t_init = pwp.times.front();
    double t_final = pwp.times.back();
    double deltaT = (t_final - t_init) / samples;

    // Initialize p_last
    geometry_msgs::msg::Point p_last = convertEigen2Point(pwp.eval(t_init));

    // Initialize j for the marker id // TODO: hardcoded
    int j = 7 * 9000;

    // Loop through the time interval
    for (double t = t_init; t <= t_final; t = t + deltaT)
    {
      visualization_msgs::msg::Marker m;
      m.type = visualization_msgs::msg::Marker::ARROW;
      m.header.frame_id = "world";
      m.header.stamp = rclcpp::Clock().now();
      m.action = visualization_msgs::msg::Marker::ADD;
      m.id = j;
      m.color = getColor(mighty_utils::red_normal);
      m.scale.x = 0.1;
      m.scale.y = 0.0000001; // rviz complains if not
      m.scale.z = 0.0000001; // rviz complains if not
      m.pose.orientation.w = 1.0;

      // Get the point
      geometry_msgs::msg::Point p = convertEigen2Point(pwp.eval(t));
      m.points.push_back(p_last);
      m.points.push_back(p);

      // Update p_last
      p_last = p;

      // Push back the marker
      marker_array.markers.push_back(m);
      j = j + 1;
    }

    // Return the marker array
    return marker_array;
  }

  geometry_msgs::msg::Point convertEigen2Point(Eigen::Vector3d vector)
  {
    geometry_msgs::msg::Point tmp;
    tmp.x = vector[0];
    tmp.y = vector[1];
    tmp.z = vector[2];
    return tmp;
  }

  // Function to convert std::vector<float> to Eigen::Vector3d
  Eigen::Vector3d convertCovMsg2Cov(const std::vector<float> &msg_cov)
  {
    // Ensure the vector has exactly 3 elements
    if (msg_cov.size() != 3)
    {
      throw std::invalid_argument("msg_cov must have exactly 3 elements.");
    }

    // Convert to Eigen::Vector3d
    Eigen::Vector3d cov;
    cov << msg_cov[0], msg_cov[1], msg_cov[2];
    return cov;
  }

  // Function to convert Eigen::Vector3d to std::vector<float>
  std::vector<float> convertCov2CovMsg(const Eigen::Vector3d &cov)
  {
    // Convert to std::vector<float>
    std::vector<float> msg_cov;
    msg_cov.push_back(cov[0]);
    msg_cov.push_back(cov[1]);
    msg_cov.push_back(cov[2]);
    return msg_cov;
  }

  // Function to convert std::vector<float> to Eigen::Matrix<double, 6, 1>
  Eigen::Matrix<double, 6, 1> convertCoeffMsg2Coeff(const std::vector<float> &msg_coeff)
  {
    // Ensure the vector has exactly 6
    if (msg_coeff.size() != 6)
    {
      throw std::invalid_argument("msg_coeff must have exactly 6 elements.");
    }

    // Convert to Eigen::Matrix<double, 6, 1>
    Eigen::Matrix<double, 6, 1> coeff;
    coeff << msg_coeff[0], msg_coeff[1], msg_coeff[2], msg_coeff[3], msg_coeff[4], msg_coeff[5]; 
    return coeff;
  }

  std_msgs::msg::ColorRGBA getColor(int id)
  {

    // Define colors

    // red
    std_msgs::msg::ColorRGBA red;
    red.r = 1;
    red.g = 0;
    red.b = 0;
    red.a = 1;

    // red transparent
    std_msgs::msg::ColorRGBA red_trans;
    red_trans.r = 1;
    red_trans.g = 0;
    red_trans.b = 0;
    red_trans.a = 0.7;

    // red transparent transparent
    std_msgs::msg::ColorRGBA red_trans_trans;
    red_trans_trans.r = 1;
    red_trans_trans.g = 0;
    red_trans_trans.b = 0;
    red_trans_trans.a = 0.4;

    // blue
    std_msgs::msg::ColorRGBA blue;
    blue.r = 0;
    blue.g = 0;
    blue.b = 1;
    blue.a = 1;

    // blue transparent
    std_msgs::msg::ColorRGBA blue_trans;
    blue_trans.r = 0;
    blue_trans.g = 0;
    blue_trans.b = 1;
    blue_trans.a = 0.7;

    // blue transparent transparent
    std_msgs::msg::ColorRGBA blue_trans_trans;
    blue_trans_trans.r = 0;
    blue_trans_trans.g = 0;
    blue_trans_trans.b = 1;
    blue_trans_trans.a = 0.4;

    // blue light
    std_msgs::msg::ColorRGBA blue_light;
    blue_light.r = 0.5;
    blue_light.g = 0.7;
    blue_light.b = 1;
    blue_light.a = 1;

    // green
    std_msgs::msg::ColorRGBA green;
    green.r = 0;
    green.g = 1;
    green.b = 0;
    green.a = 1;

    // yellow
    std_msgs::msg::ColorRGBA yellow;
    yellow.r = 1;
    yellow.g = 1;
    yellow.b = 0;
    yellow.a = 1;

    // orange transparent
    std_msgs::msg::ColorRGBA orange_trans; // orange transparent
    orange_trans.r = 1;
    orange_trans.g = 0.5;
    orange_trans.b = 0;
    orange_trans.a = 0.7;

    // teal normal
    std_msgs::msg::ColorRGBA teal_normal; // orange transparent
    teal_normal.r = 25 / 255.0;
    teal_normal.g = 1.0;
    teal_normal.b = 240.0 / 255.0;
    teal_normal.a = 1.0;

    // black transparent
    std_msgs::msg::ColorRGBA black_trans; // orange transparent
    black_trans.r = 0.0;
    black_trans.g = 0.0;
    black_trans.b = 0.0;
    black_trans.a = 0.7;

    std_msgs::msg::ColorRGBA green_trans_trans;  // green transparent
    green_trans_trans.r = 0;
    green_trans_trans.g = 1;
    green_trans_trans.b = 0;
    green_trans_trans.a = 0.2;

    switch (id)
    {
    case mighty_utils::red_normal:
      return red;
      break;
    case mighty_utils::red_trans:
      return red_trans;
      break;
    case mighty_utils::red_trans_trans:
      return red_trans_trans;
      break;
    case mighty_utils::blue_normal:
      return blue;
      break;
    case mighty_utils::blue_trans:
      return blue_trans;
      break;
    case mighty_utils::blue_trans_trans:
      return blue_trans_trans;
      break;
    case mighty_utils::blue_light:
      return blue_light;
      break;
    case mighty_utils::green_normal:
      return green;
      break;
    case mighty_utils::yellow_normal:
      return yellow;
      break;
    case mighty_utils::orange_trans:
      return orange_trans;
      break;
    case mighty_utils::black_trans:
      return black_trans;
      break;
    case mighty_utils::teal_normal:
      return teal_normal;
      break;
    case mighty_utils::green_trans_trans:
      return green_trans_trans;
      break;
    default:
      std::cout << "COLOR NOT DEFINED, returning RED" << std::endl;
      return red;
    }
  }

  std::vector<Eigen::Matrix<double, 3, 4>> convertCoefficients2ControlPoints(const PieceWisePol &pwp, const Eigen::Matrix<double, 4, 4> &A_rest_pos_basis_inverse)
  {

    // Initialization of the control points
    std::vector<Eigen::Matrix<double, 3, 4>> control_points;

    for (int i = 0; i < pwp.coeff_x.size(); i++)
    {
      // Get P matrix
      Eigen::Matrix<double, 3, 4> P;
      P.row(0) = pwp.coeff_x[i];
      P.row(1) = pwp.coeff_y[i];
      P.row(2) = pwp.coeff_z[i];

      // Get V matrix
      Eigen::Matrix<double, 3, 4> V = P * A_rest_pos_basis_inverse;

      // Push back the control points
      control_points.push_back(V);
    }

    return control_points;
  }

  double getMinTimeDoubleIntegrator1D(const double p0, const double v0, const double pf, const double vf,
                                      const double v_max, const double a_max)
  {
    // The notation of this function is based on the paper "Constrained time-optimal control of double integrator system and its application in MPC"
    // https://iopscience.iop.org/article/10.1088/1742-6596/783/1/012024

    double x1 = v0;
    double x2 = p0;
    double x1r = vf;
    double x2r = pf;

    double k1 = a_max; // Note that the paper uses u\in[-1, 1].But setting k1 to a_max has the same effect k2 = 1.0;
    double k2 = 1.0;

    double x1_bar = v_max;

    double B = (k2 / (2 * k1)) * sgn(-x1 + x1r) * (pow(x1, 2) - pow(x1r, 2)) + x2r;
    double C = (k2 / (2 * k1)) * (pow(x1, 2) + pow(x1r, 2)) - (k2 / k1) * pow(x1_bar, 2) + x2r;
    double D = (-k2 / (2 * k1)) * (pow(x1, 2) + pow(x1r, 2)) + (k2 / k1) * pow(x1_bar, 2) + x2r;

    double time;

    if ((x2 <= B) && (x2 >= C))
    {
      time = (-k2 * (x1 + x1r) +
              2 * sqrt(pow(k2, 2) * pow(x1, 2) - k1 * k2 * ((k2 / (2 * k1)) * (pow(x1, 2) - pow(x1r, 2)) + x2 - x2r))) /
             (k1 * k2);
    }

    else if ((x2 <= B) && (x2 < C))
    {
      time = (x1_bar - x1 - x1r) / k1 + (pow(x1, 2) + pow(x1r, 2)) / (2 * k1 * x1_bar) + (x2r - x2) / (k2 * x1_bar);
    }

    else if ((x2 > B) && (x2 <= D))
    {
      time = (k2 * (x1 + x1r) +
              2 * sqrt(pow(k2, 2) * pow(x1, 2) + k1 * k2 * ((k2 / (2 * k1)) * (-pow(x1, 2) + pow(x1r, 2)) + x2 - x2r))) /
             (k1 * k2);
    }

    else
    { // (x2 > B) && (x2 > D)

      time = (x1_bar + x1 + x1r) / k1 + (pow(x1, 2) + pow(x1r, 2)) / (2 * k1 * x1_bar) + (-x2r + x2) / (k2 * x1_bar);
    }

    return time;
  }

  double getMinTimeDoubleIntegrator3D(const Eigen::Vector3d &p0, const Eigen::Vector3d &v0, const Eigen::Vector3d &pf,
                                      const Eigen::Vector3d &vf, const Eigen::Vector3d &v_max,
                                      const Eigen::Vector3d &a_max)
  {
    double min_x = getMinTimeDoubleIntegrator1D(p0.x(), v0.x(), pf.x(), vf.x(), v_max.x(), a_max.x());
    double min_y = getMinTimeDoubleIntegrator1D(p0.y(), v0.y(), pf.y(), vf.y(), v_max.y(), a_max.y());
    double min_z = getMinTimeDoubleIntegrator1D(p0.z(), v0.z(), pf.z(), vf.z(), v_max.z(), a_max.z());

    double min_time = std::max({min_x, min_y, min_z}); // Note that it's the maximum of all the axes

    return min_time;
  }

  void quaternion2Euler(tf2::Quaternion q, double &roll, double &pitch, double &yaw)
  {
    tf2::Matrix3x3(q).getRPY(roll, pitch, yaw);
  }

  void quaternion2Euler(Eigen::Quaterniond q, double &roll, double &pitch, double &yaw)
  {
    tf2::Quaternion tf_q(q.x(), q.y(), q.z(), q.w());
    quaternion2Euler(tf_q, roll, pitch, yaw);
  }

  void quaternion2Euler(geometry_msgs::msg::Quaternion q, double &roll, double &pitch, double &yaw)
  {
    tf2::Quaternion tf_q(q.x, q.y, q.z, q.w);
    quaternion2Euler(tf_q, roll, pitch, yaw);
  }

  void saturate(double &var, double min, double max)
  {
    if (var < min)
    {
      var = min;
    }
    else if (var > max)
    {
      var = max;
    }
  }

  void angle_wrap(double &diff)
  {
    diff = fmod(diff + M_PI, 2 * M_PI);
    if (diff < 0)
      diff += 2 * M_PI;
    diff -= M_PI;
  }

  // P1-P2 is the direction used for projection. P2 is the terminal goal. wdx, wdy and wdz are the widths of a 3D box
  // centered on P1
  Eigen::Vector3d projectPointToBox(Eigen::Vector3d &P1, Eigen::Vector3d &P2, double wdx, double wdy, double wdz)
  {
    double x_max = P1(0) + wdx / 2;
    double x_min = P1(0) - wdx / 2;
    double y_max = P1(1) + wdy / 2;
    double y_min = P1(1) - wdy / 2;
    double z_max = P1(2) + wdz / 2;
    double z_min = P1(2) - wdz / 2;

    if ((P2(0) < x_max && P2(0) > x_min) && (P2(1) < y_max && P2(1) > y_min) && (P2(2) < z_max && P2(2) > z_min))
    {
      // goal is inside the map
      return P2;
    }

    // goal is outside the map
    Eigen::Vector3d inters;
    std::vector<Eigen::Vector4d> all_planes = {
        Eigen::Vector4d(1, 0, 0, -x_max), // Plane X right
        Eigen::Vector4d(-1, 0, 0, x_min), // Plane X left
        Eigen::Vector4d(0, 1, 0, -y_max), // Plane Y right
        Eigen::Vector4d(0, -1, 0, y_min), // Plane Y left
        Eigen::Vector4d(0, 0, 1, -z_max), // Plane Z up
        Eigen::Vector4d(0, 0, -1, z_min)  // Plane Z down
    };

    vec_Vecf<3> intersections;

    int axis; // 1 is x, 2 is y, 3 is z
    for (int i = 0; i < 6; i++)
    {
      if (getIntersectionWithPlane(P1, P2, all_planes[i], inters) == true)
      {
        intersections.push_back(inters);
      }
    }

    if (intersections.size() == 0)
    {
      // There is no intersection
      RCLCPP_ERROR(rclcpp::get_logger("rclcpp"), "This is impossible, there should be an intersection");
    }

    // Initialize the distances
    std::vector<double> distances;
    for (size_t i = 0; i < intersections.size(); i++)
    {
      double distance = (intersections[i] - P1).norm();
      distances.push_back(distance);
    }

    // Get the minimum element index
    int minElementIndex = std::min_element(distances.begin(), distances.end()) - distances.begin();
    inters = intersections[minElementIndex];

    // Sometimes observe that the goal is projected just outside the box. This is a hack to make sure the goal is inside
    // Move inters to the direction (from inters to P1) by a small amount
    Eigen::Vector3d v = (inters - P1).normalized();
    inters = inters - v * 1.0;

    return inters;
  }

  // Project point to a sphere
  Eigen::Vector3d projectPointToSphere(const Eigen::Vector3d &P1, const Eigen::Vector3d &P2, double radius)
  {

    // If the point is already inside the sphere, return the point
    if ((P2 - P1).norm() <= radius)
    {
      return P2;
    }

    Eigen::Vector3d inters;
    Eigen::Vector3d v = (P2 - P1).normalized();
    inters = P1 + v * radius;
    return inters;
  }

  // Convert visualization_msgs::msg::MarkerArray to vec_Vecf<3>
  void convertMarkerArray2Vec_Vec_Vecf3(const visualization_msgs::msg::MarkerArray &marker_array, std::vector<vec_Vecf<3>> &vec, std::vector<double> &scale_vec)
  {
    // loop through the markers
    for (int i = 0; i < marker_array.markers.size(); i++)
    {
      // if the marker is empty, skip
      if (marker_array.markers[i].points.size() == 0)
      {
        continue;
      }

      // push back the scale
      scale_vec.push_back(marker_array.markers[i].scale.x);

      // initialize the vector of points
      vec_Vecf<3> vec_tmp;

      // loop through the points
      for (int j = 0; j < marker_array.markers[i].points.size(); j++)
      {
        Eigen::Vector3d point;
        point << marker_array.markers[i].points[j].x, marker_array.markers[i].points[j].y, marker_array.markers[i].points[j].z;
        vec_tmp.push_back(point);
      }
      vec.push_back(vec_tmp);
    }
  }

  void createMoreVertexes(vec_Vecf<3> &path, double d)
  {

    if (path.size() < 2)
    {
      return;
    }

    for (int j = 0; j < path.size() - 1; j++)
    {
      double dist = (path[j + 1] - path[j]).norm();
      int vertexes_to_add = floor(dist / d);
      Eigen::Vector3d v = (path[j + 1] - path[j]).normalized();
      if (dist > d)
      {
        for (int i = 0; i < vertexes_to_add; i++)
        {
          path.insert(path.begin() + j + 1, path[j] + v * d);
          j = j + 1;
        }
      }
    }
  }

  double euclideanDistance(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2)
  {
    return std::sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
                     (p1[1] - p2[1]) * (p1[1] - p2[1]) +
                     (p1[2] - p2[2]) * (p1[2] - p2[2]));
  }

  Eigen::Matrix4d transformStampedToMatrix(const geometry_msgs::msg::TransformStamped &transform_stamped)
  {
    // Extract translation
    double tx = transform_stamped.transform.translation.x;
    double ty = transform_stamped.transform.translation.y;
    double tz = transform_stamped.transform.translation.z;

    // Extract rotation (quaternion)
    double qx = transform_stamped.transform.rotation.x;
    double qy = transform_stamped.transform.rotation.y;
    double qz = transform_stamped.transform.rotation.z;
    double qw = transform_stamped.transform.rotation.w;

    // Convert quaternion to rotation matrix
    Eigen::Quaterniond quaternion(qw, qx, qy, qz);
    Eigen::Matrix3d rotation_matrix = quaternion.toRotationMatrix();

    // Build the 4x4 transformation matrix
    Eigen::Matrix4d transformation_matrix = Eigen::Matrix4d::Identity();
    transformation_matrix.block<3, 3>(0, 0) = rotation_matrix; // Top-left 3x3 rotation
    transformation_matrix(0, 3) = tx;                          // Top-right x translation
    transformation_matrix(1, 3) = ty;                          // Top-right y translation
    transformation_matrix(2, 3) = tz;                          // Top-right z translation

    return transformation_matrix;
  }

  /**
   * @brief Find velocity at each point in the path
   * @param path vec_Vecf<3>
   * @param velocities vec_Vecf<3>
   * @param A state
   */
  void findVelocitiesInPath(const vec_Vecf<3> &path, vec_Vecf<3> &velocities, const state &A, const Eigen::Vector3d &v_max_3d, bool verbose)
  {

    // Get the initial velocity (which is A.vel)
    velocities.push_back(A.vel);

    // Threshold distance for max velocity saturation
    double threshold_distance = 5.0; // TODO: make this a parameter

    // Loop through the global path to be used in
    for (int i = 0; i < path.size() - 2; i++)
    {

      // Get point 0, 1, and 2
      Eigen::Vector3d point0(path[i].x(), path[i].y(), path[i].z());
      Eigen::Vector3d point1(path[i + 1].x(), path[i + 1].y(), path[i + 1].z());
      Eigen::Vector3d point2(path[i + 2].x(), path[i + 2].y(), path[i + 2].z());

      // Compute direction from point 0 to point 1
      Eigen::Vector3d direction01 = point1 - point0;

      // Compute direction from point 1 to point 2
      Eigen::Vector3d direction12 = point2 - point1;

      // Compute direction from point 0 to point 2
      Eigen::Vector3d direction02 = point2 - point0;

      // Initialize velocity_{i+1}
      Eigen::Vector3d velocity_i_plus_1;

      // Loop through x, y, and z
      for (int j = 0; j < 3; j++)
      {

        // Case 1: positive + positive
        if (direction01[j] > 0 && direction12[j] > 0 || direction01[j] == 0 && direction12[j] > 0 || direction01[j] > 0 && direction12[j] == 0)
        {
          // Compute the velocity based on the distance between point 0 and point 2
          double dist02 = abs(direction02[j]); // distance on axis j between point 0 and point 2
          if (dist02 > threshold_distance)
          {
            // If the distance is bigger than the threshold distance, use the max velocity
            velocity_i_plus_1[j] = v_max_3d[j];
          }
          else
          {
            // If the distance is smaller than the threshold distance, use the velocity based on the distance
            velocity_i_plus_1[j] = v_max_3d[j] * dist02 / threshold_distance;
          }

          continue;
        } // End of Case 1

        // Case 2: positive and negative
        if (direction01[j] > 0 && direction12[j] < 0)
        {
          // The agent needs to change direction -> set the velocity to 0
          velocity_i_plus_1[j] = 0.0;

          continue;
        } // End of Case 2

        // Case 3: negative and positive
        if (direction01[j] < 0 && direction12[j] > 0)
        {
          // The agent needs to change direction -> set the velocity to 0
          velocity_i_plus_1[j] = 0.0;

          continue;
        } // End of Case 3

        // Case 4: negative and negative
        if (direction01[j] < 0 && direction12[j] < 0 || direction01[j] == 0 && direction12[j] < 0 || direction01[j] < 0 && direction12[j] == 0)
        {
          // Compute the velocity based on the distance between point 0 and point 2
          double dist02 = abs(direction02[j]); // distance on axis j between point 0 and point 2
          if (dist02 > threshold_distance)
          {
            // If the distance is bigger than the threshold distance, use the max velocity
            velocity_i_plus_1[j] = -v_max_3d[j];
          }
          else
          {
            // If the distance is smaller than the threshold distance, use the velocity based on the distance
            velocity_i_plus_1[j] = -v_max_3d[j] * dist02 / threshold_distance;
          }

          continue;
        } // End of Case 4

        // Case 5: zero ans zero
        if (direction01[j] == 0 && direction12[j] == 0)
        {
          // The agent is not moving -> set the velocity to 0
          velocity_i_plus_1[j] = 0.0;

          continue;
        } // End of Case 5

      } // End of loop through x, y, and z

      // Add the velocity to the vector
      velocities.push_back(velocity_i_plus_1);

    } // End of loop through the global path to be used in

    // Add the last velocity (which is zero)
    velocities.push_back(Eigen::Vector3d::Zero());

    // Debug
    assert(path.size() == velocities.size() && "path.size() != velocities.size()");

    // print out the velocities
    // if (verbose)
    // {
    //   for (int i = 0; i < velocities.size(); i++)
    //   {
    //     std::cout << "Velocity " << i << ": " << velocities[i].transpose() << std::endl;
    //   }
    // }

  } // End of findVelocitiesInPath

  /**
   * @brief Get time allocation for each polytope
   * @param path vec_Vecf<3>
   * @param A state
   * @return travel_times std::vector<double>
   */
  std::vector<double> getTravelTimes(const vec_Vecf<3> &path, const state &A, bool debug_verbose, const Eigen::Vector3d &v_max_3d, const Eigen::Vector3d &a_max_3d)
  {

    // First, estimate velocity at each point
    vec_Vecf<3> velocities;
    velocities.reserve(path.size());
    if (path.size() == 1)
    {
      return {};
    }
    else if (path.size() == 2)
    {
      double travel_time = mighty_utils::getMinTimeDoubleIntegrator3D(path[0], A.vel, path[1], Eigen::Vector3d::Zero(), v_max_3d, a_max_3d);
      return {travel_time};
    }
    else
    {
      findVelocitiesInPath(path, velocities, A, v_max_3d, debug_verbose);
    }

    // Second, find the minimum time allocation

    // Initialize travel times
    std::vector<double> travel_times;

    // Loop through the initial guess path and find travel times
    for (int i = 0; i < path.size() - 1; i++)
    {

      // Compute travel time
      double travel_time = mighty_utils::getMinTimeDoubleIntegrator3D(path[i], velocities[i], path[i + 1], velocities[i + 1], v_max_3d, a_max_3d);

      if (debug_verbose)
      {
        std::cout << "i: " << i << " of " << path.size() - 1 << std::endl;
        std::cout << "Start: " << path[i].transpose() << std::endl;
        std::cout << "start velocity: " << velocities[i].transpose() << std::endl;
        std::cout << "End: " << path[i + 1].transpose() << std::endl;
        std::cout << "end velocity: " << velocities[i + 1].transpose() << std::endl;
        std::cout << "travel time: " << travel_time << std::endl;
      }
      // Add travel time to the vector
      travel_times.push_back(travel_time);
    }

    return travel_times;

  } // End of getTravelTimes

  geometry_msgs::msg::Pose identityGeometryMsgsPose()
  {
    geometry_msgs::msg::Pose pose;
    pose.position.x = 0;
    pose.position.y = 0;
    pose.position.z = 0;
    pose.orientation.x = 0;
    pose.orientation.y = 0;
    pose.orientation.z = 0;
    pose.orientation.w = 1;
    return pose;
  }

} // namespace mighty_utils