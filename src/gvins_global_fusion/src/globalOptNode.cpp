/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 *
 * This file is part of VINS.
 *
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Qin Tong (qintonguav@gmail.com)
 *******************************************************/

#include "globalOpt.h"
#include <atomic>
#include <gnss_comm/GnssPVTSolnMsg.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <queue>
#include <ros/ros.h>
#include <sensor_msgs/NavSatFix.h>
#include <stdio.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

GlobalOptimization globalEstimator;
double GNSS_LOCAL_TIME_DIFF = 18.0;
ros::Publisher pub_global_odometry, pub_global_path, pub_car;
nav_msgs::Path *global_path;
double last_vio_t = -1;
std::queue<sensor_msgs::NavSatFixConstPtr> gpsQueue;
std::mutex m_buf;

static bool g_verbose = false;
static bool g_save_csv = false;
static std::string g_mesh_package = "gvins_global_fusion";
static std::string g_csv_path;
static double g_sync_tolerance = 0.05;
static std::ofstream g_csv_stream;
static bool g_csv_header_written = false;
/** 255 = unknown; carr_soln 0/1/2 per GnssPVTSolnMsg (2 = RTK carrier-phase fix). */
static std::atomic<uint8_t> g_carr_soln{255};
static bool g_use_pvt_gate = false;

static void pvt_callback(const gnss_comm::GnssPVTSolnMsgConstPtr &msg)
{
    g_carr_soln.store(msg->carr_soln);
}

void publish_car_model(double t, Eigen::Vector3d t_w_car, Eigen::Quaterniond q_w_car)
{
    visualization_msgs::MarkerArray markerArray_msg;
    visualization_msgs::Marker car_mesh;
    car_mesh.header.stamp = ros::Time(t);
    car_mesh.header.frame_id = "world";
    car_mesh.type = visualization_msgs::Marker::MESH_RESOURCE;
    car_mesh.action = visualization_msgs::Marker::ADD;
    car_mesh.id = 0;

    car_mesh.mesh_resource = "package://" + g_mesh_package + "/models/car.dae";

    Eigen::Matrix3d rot;
    rot << 0, 0, -1, 0, -1, 0, -1, 0, 0;

    Eigen::Quaterniond Q;
    Q = q_w_car * rot;
    car_mesh.pose.position.x = t_w_car.x();
    car_mesh.pose.position.y = t_w_car.y();
    car_mesh.pose.position.z = t_w_car.z();
    car_mesh.pose.orientation.w = Q.w();
    car_mesh.pose.orientation.x = Q.x();
    car_mesh.pose.orientation.y = Q.y();
    car_mesh.pose.orientation.z = Q.z();

    car_mesh.color.a = 1.0;
    car_mesh.color.r = 1.0;
    car_mesh.color.g = 0.0;
    car_mesh.color.b = 0.0;

    float major_scale = 2.0;

    car_mesh.scale.x = major_scale;
    car_mesh.scale.y = major_scale;
    car_mesh.scale.z = major_scale;
    markerArray_msg.markers.push_back(car_mesh);
    pub_car.publish(markerArray_msg);
}

void GPS_callback(const sensor_msgs::NavSatFixConstPtr &GPS_msg)
{
    m_buf.lock();
    gpsQueue.push(GPS_msg);
    m_buf.unlock();
}

void vio_callback(const nav_msgs::Odometry::ConstPtr &pose_msg)
{
    double t = pose_msg->header.stamp.toSec();
    last_vio_t = t;
    Eigen::Vector3d vio_t(pose_msg->pose.pose.position.x, pose_msg->pose.pose.position.y,
                          pose_msg->pose.pose.position.z);
    Eigen::Quaterniond vio_q;
    vio_q.w() = pose_msg->pose.pose.orientation.w;
    vio_q.x() = pose_msg->pose.pose.orientation.x;
    vio_q.y() = pose_msg->pose.pose.orientation.y;
    vio_q.z() = pose_msg->pose.pose.orientation.z;
    globalEstimator.inputOdom(t, vio_t, vio_q);

    m_buf.lock();
    double t_aligned = t + GNSS_LOCAL_TIME_DIFF;
    while (!gpsQueue.empty())
    {
        sensor_msgs::NavSatFixConstPtr GPS_msg = gpsQueue.front();
        double gps_t = GPS_msg->header.stamp.toSec();
        if (g_verbose)
            printf("vio t: %f, vio t aligned: %f, gps t: %f \n", t, t_aligned, gps_t);
        if (gps_t >= t_aligned - g_sync_tolerance && gps_t <= t_aligned + g_sync_tolerance)
        {
            double latitude = GPS_msg->latitude;
            double longitude = GPS_msg->longitude;
            double altitude = GPS_msg->altitude;
            double pos_accuracy = GPS_msg->position_covariance[0];
            if (pos_accuracy <= 0)
                pos_accuracy = 1;
            const bool loose_ok = (!g_use_pvt_gate || g_carr_soln.load() == 2);
            if (loose_ok)
                globalEstimator.inputGPS(t, latitude, longitude, altitude, pos_accuracy);
            gpsQueue.pop();
            break;
        }
        else if (gps_t < t_aligned - g_sync_tolerance)
            gpsQueue.pop();
        else if (gps_t > t_aligned + g_sync_tolerance)
            break;
    }
    m_buf.unlock();

    Eigen::Vector3d global_t;
    Eigen::Quaterniond global_q;
    globalEstimator.getGlobalOdom(t, global_t, global_q);
    global_q.normalize();

    nav_msgs::Odometry odometry;
    odometry.header = pose_msg->header;
    odometry.header.frame_id = "world";
    odometry.child_frame_id = "world";
    odometry.pose.pose.position.x = global_t.x();
    odometry.pose.pose.position.y = global_t.y();
    odometry.pose.pose.position.z = global_t.z();
    odometry.pose.pose.orientation.x = global_q.x();
    odometry.pose.pose.orientation.y = global_q.y();
    odometry.pose.pose.orientation.z = global_q.z();
    odometry.pose.pose.orientation.w = global_q.w();
    pub_global_odometry.publish(odometry);
    pub_global_path.publish(*global_path);
    publish_car_model(t, global_t, global_q);

    if (g_save_csv && g_csv_stream.is_open())
    {
        if (!g_csv_header_written)
        {
            g_csv_stream << "# timestamp_s x y z roll pitch yaw yaw_rate\n";
            g_csv_header_written = true;
        }
        g_csv_stream.setf(std::ios::fixed, std::ios::floatfield);
        g_csv_stream.precision(16);
        g_csv_stream << pose_msg->header.stamp.toSec() + GNSS_LOCAL_TIME_DIFF << " ";
        g_csv_stream.precision(5);
        g_csv_stream << global_t.x() << " " << global_t.y() << " " << global_t.z() << " "
                      << "0.0"
                      << " "
                      << "0.0"
                      << " "
                      << "0.0"
                      << " "
                      << "0.0" << std::endl;
    }
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "gvins_global_fusion");
    ros::NodeHandle nh("~");

    global_path = &globalEstimator.global_path;

    std::string vio_odom_topic;
    std::string gps_topic;
    nh.param<std::string>("vio_odom_topic", vio_odom_topic, std::string("/gvins/odometry"));
    nh.param<std::string>("gps_topic", gps_topic, std::string("/gps"));
    nh.param("gnss_local_time_diff", GNSS_LOCAL_TIME_DIFF, 18.0);
    nh.param("sync_tolerance", g_sync_tolerance, 0.05);
    nh.param("verbose", g_verbose, false);
    nh.param("save_csv", g_save_csv, false);
    nh.param<std::string>("csv_path", g_csv_path, std::string(""));
    nh.param<std::string>("mesh_package", g_mesh_package, std::string("gvins_global_fusion"));
    nh.param("use_pvt_gate", g_use_pvt_gate, false);
    std::string pvt_topic;
    nh.param<std::string>("pvt_topic", pvt_topic, std::string(""));
    ros::Subscriber sub_pvt;
    if (g_use_pvt_gate && !pvt_topic.empty())
    {
        sub_pvt = nh.subscribe(pvt_topic, 500, pvt_callback, ros::TransportHints().tcpNoDelay());
        ROS_INFO("gvins_global_fusion: PVT gate ON - NavSatFix fused only when carr_soln==2 (RTK carrier fix). topic=%s",
                 pvt_topic.c_str());
    }

    if (g_save_csv)
    {
        if (g_csv_path.empty())
        {
            ROS_WARN("save_csv is true but csv_path is empty; CSV logging disabled.");
            g_save_csv = false;
        }
        else
        {
            g_csv_stream.open(g_csv_path.c_str(), std::ios::out | std::ios::trunc);
            if (!g_csv_stream.is_open())
            {
                ROS_WARN("Could not open csv_path '%s'; CSV logging disabled.", g_csv_path.c_str());
                g_save_csv = false;
            }
        }
    }

    ROS_INFO("gvins_global_fusion: vio_odom_topic=%s gps_topic=%s gnss_local_time_diff=%.3f sync_tolerance=%.3f",
             vio_odom_topic.c_str(), gps_topic.c_str(), GNSS_LOCAL_TIME_DIFF, g_sync_tolerance);

    // Large queues: VIO is often ~10–20 Hz; 100 messages is only ~5–10 s of backlog — under load
    // or brief stalls, drops here truncate global_odometry / CSV / hybrid TUM on long bags.
    ros::Subscriber sub_GPS =
        nh.subscribe(gps_topic, 2000, GPS_callback, ros::TransportHints().tcpNoDelay());
    ros::Subscriber sub_vio =
        nh.subscribe(vio_odom_topic, 2000, vio_callback, ros::TransportHints().tcpNoDelay());

    pub_global_path = nh.advertise<nav_msgs::Path>("global_path", 1000);
    pub_global_odometry = nh.advertise<nav_msgs::Odometry>("global_odometry", 1000);
    pub_car = nh.advertise<visualization_msgs::MarkerArray>("car_model", 1000);
    ros::spin();
    return 0;
}
