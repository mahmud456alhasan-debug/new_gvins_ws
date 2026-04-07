/**
 * Single-frame TUM for evo vs RTK GT: logs gvins_global_fusion/global_odometry only.
 *
 * Adaptive tight/loose still changes GVINS before this node; global fusion consumes /gvins/odometry,
 * so the exported trajectory stays in one GPS-aligned frame (same idea as loose_global_fusion.csv).
 *
 * Timestamps: stamp + gnss_local_time_diff (match gvins_global_fusion loose CSV).
 *
 * Legacy mode log_mode:=split mixed /gvins/odometry with global odometry — different frames; do not use
 * with evo --align_origin (use evo_ape -a or avoid split).
 */

#include <fstream>
#include <iomanip>
#include <string>

#include <geometry_msgs/Pose.h>
#include <nav_msgs/Odometry.h>
#include <ros/ros.h>

namespace
{
std::ofstream *g_out{nullptr};
double g_time_diff{0.0};

void write_tum_line(double t, const geometry_msgs::Pose &p)
{
    if (!g_out || !g_out->good())
        return;
    *g_out << std::setprecision(18) << t << " " << p.position.x << " " << p.position.y << " " << p.position.z << " "
           << p.orientation.x << " " << p.orientation.y << " " << p.orientation.z << " " << p.orientation.w << "\n";
    g_out->flush();
}

void global_only_cb(const nav_msgs::Odometry::ConstPtr &msg)
{
    const double t = msg->header.stamp.toSec() + g_time_diff;
    write_tum_line(t, msg->pose.pose);
}
} // namespace

int main(int argc, char **argv)
{
    ros::init(argc, argv, "hybrid_trajectory_logger");
    ros::NodeHandle nh("~");

    std::string out_path;
    std::string global_odom_topic;
    std::string log_mode;

    nh.param<std::string>("output_path", out_path, std::string(""));
    nh.param<std::string>("global_odom_topic", global_odom_topic, std::string("/gvins_global_fusion/global_odometry"));
    nh.param("gnss_local_time_diff", g_time_diff, 18.0);
    nh.param<std::string>("log_mode", log_mode, std::string("global_only"));

    if (out_path.empty())
    {
        ROS_ERROR("Set ~output_path, e.g. $HOME/new_gvins_ws/output/hybrid_trajectory_tum.txt");
        return 1;
    }

    std::ofstream out(out_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.is_open())
    {
        ROS_ERROR("Cannot open %s", out_path.c_str());
        return 1;
    }
    g_out = &out;

    if (log_mode == "split")
    {
        ROS_ERROR("log_mode=split mixes coordinate frames; not supported in this build. Use global_only.");
        return 1;
    }

    ros::Subscriber sub_glob =
        nh.subscribe<nav_msgs::Odometry>(global_odom_topic, 2000, global_only_cb, ros::TransportHints().tcpNoDelay());

    ROS_INFO("Writing TUM to %s (single frame: %s, t += gnss_local_time_diff=%.3f s)", out_path.c_str(),
             global_odom_topic.c_str(), g_time_diff);

    ros::spin();
    g_out = nullptr;
    out.close();
    return 0;
}
