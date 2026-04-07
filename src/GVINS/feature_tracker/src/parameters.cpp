#include "parameters.h"

std::string IMAGE_TOPIC;
std::string IMU_TOPIC;
std::vector<std::string> CAM_NAMES;
std::string FISHEYE_MASK;
int MAX_CNT;
int MIN_DIST;
int WINDOW_SIZE;
int FREQ;
double F_THRESHOLD;
int SHOW_TRACK;
int STEREO_TRACK;
int EQUALIZE;
int ROW;
int COL;
int FOCAL_LENGTH;
int FISHEYE;
bool PUB_THIS_FRAME;

template <typename T>
T readParam(ros::NodeHandle &n, std::string name)
{
    T ans;
    if (n.getParam(name, ans))
    {
        ROS_INFO_STREAM("Loaded " << name << ": " << ans);
    }
    else
    {
        ROS_ERROR_STREAM("Failed to load " << name);
        n.shutdown();
    }
    return ans;
}

void readParameters(ros::NodeHandle &n)
{
    std::string config_file;
    config_file = readParam<std::string>(n, "config_file");
    cv::FileStorage fsSettings(config_file, cv::FileStorage::READ);
    if(!fsSettings.isOpened())
    {
        std::cerr << "ERROR: Wrong path to settings" << std::endl;
    }
    std::string GVINS_FOLDER_PATH = readParam<std::string>(n, "gvins_folder");

    fsSettings["image_topic"] >> IMAGE_TOPIC;
    fsSettings["imu_topic"] >> IMU_TOPIC;
    MAX_CNT = fsSettings["max_cnt"];
    MIN_DIST = fsSettings["min_dist"];
    ROW = fsSettings["image_height"];
    COL = fsSettings["image_width"];
    FREQ = fsSettings["freq"];
    F_THRESHOLD = fsSettings["F_threshold"];
    SHOW_TRACK = fsSettings["show_track"];
    EQUALIZE = fsSettings["equalize"];
    FISHEYE = fsSettings["fisheye"];
    if (FISHEYE == 1)
        FISHEYE_MASK = GVINS_FOLDER_PATH + "/config/fisheye_mask.jpg";
    CAM_NAMES.push_back(config_file);

    WINDOW_SIZE = 20;
    STEREO_TRACK = false;
    FOCAL_LENGTH = 460;
    PUB_THIS_FRAME = false;

    if (FREQ == 0)
        FREQ = 100;

    // Adaptive parameter selection: Use optimized values for VIO-only mode
    // When GNSS is disabled, override to optimized parameters for better VIO performance
    int gnss_enable_value = fsSettings["gnss_enable"];
    if (gnss_enable_value == 0)  // VIO-only mode
    {
        MAX_CNT = 200;  // Optimized for VIO-only (better tracking in challenging scenes)
        MIN_DIST = 25;  // Optimized for VIO-only (more distributed features)
        ROS_INFO_STREAM("GNSS disabled: Using optimized VIO parameters (max_cnt=200, min_dist=25)");
    }
    else  // GNSS-enabled mode
    {
        // Use config values as-is (baseline values that work with GNSS initialization)
        ROS_INFO_STREAM("GNSS enabled: Using baseline parameters from config (max_cnt=" << MAX_CNT << ", min_dist=" << MIN_DIST << ")");
    }

    fsSettings.release();


}
