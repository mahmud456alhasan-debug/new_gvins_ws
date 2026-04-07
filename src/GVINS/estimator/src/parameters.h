/*
 * @Author: Demo wxtcon@163.com
 * @Date: 2023-10-09 21:58:59
 * @LastEditors: Demo wxtcon@163.com
 * @LastEditTime: 2023-11-27 19:01:50
 * @FilePath: /GVINS/estimator/src/parameters.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
 #pragma once
 
 #include <atomic>
 #include <cstdlib>
 #include <ros/ros.h>
 #include <vector>
 #include <eigen3/Eigen/Dense>
 #include "utility/utility.h"
 #include <opencv2/opencv.hpp>
 #include <opencv2/core/eigen.hpp>
 #include <fstream>
 
 const double FOCAL_LENGTH = 460.0; //焦距
 const int WINDOW_SIZE = 10;
 const int NUM_OF_CAM = 1;
 const int NUM_OF_F = 1000;
 //#define UNIT_SPHERE_ERROR
 
 extern double INIT_DEPTH;
 extern double MIN_PARALLAX;
 extern int ESTIMATE_EXTRINSIC;
 
 extern double ACC_N, ACC_W;
 extern double GYR_N, GYR_W;
 
 extern std::vector<Eigen::Matrix3d> RIC;
 extern std::vector<Eigen::Vector3d> TIC;
 extern Eigen::Vector3d G;
 
 extern double BIAS_ACC_THRESHOLD;
 extern double BIAS_GYR_THRESHOLD;
 extern double SOLVER_TIME;
 extern int NUM_ITERATIONS;
 extern std::string EX_CALIB_RESULT_PATH;
 extern std::string VINS_RESULT_PATH;
 extern std::string FACTOR_GRAPH_RESULT_PATH;
 // record file relate (damon)
 extern std::string RECORD_FILE_PATH;
 
 
 extern std::string IMU_TOPIC;
 extern double TD;
 extern int ESTIMATE_TD;
 extern double ROW, COL;
 
 extern std::atomic<bool> GNSS_ENABLE;
 /** When true: subscribe to gnss_pvt_topic and set GNSS_ENABLE at runtime — RTK fix (carr_soln==2) => loose only; carr_soln<2 => tight raw GNSS. */
 extern bool GNSS_ADAPTIVE_CARR_SOLN_SWITCH;
 extern std::string GNSS_PVT_TOPIC;
 extern std::string GNSS_EPHEM_TOPIC;
 extern std::string GNSS_GLO_EPHEM_TOPIC;
 extern std::string GNSS_MEAS_TOPIC;
 extern std::string GNSS_IONO_PARAMS_TOPIC;
 extern std::string GNSS_TP_INFO_TOPIC;
 extern std::vector<double> GNSS_IONO_DEFAULT_PARAMS;
 extern bool GNSS_LOCAL_ONLINE_SYNC;
 extern std::string LOCAL_TRIGGER_INFO_TOPIC;
 extern double GNSS_LOCAL_TIME_DIFF;
 extern double GNSS_ELEVATION_THRES;
 extern double GNSS_PSR_STD_THRES;
 extern double GNSS_DOPP_STD_THRES;
 extern uint32_t GNSS_TRACK_NUM_THRES;
 extern double GNSS_DDT_WEIGHT;
 extern std::string GNSS_RESULT_PATH;
 
 extern double GNSS_DENIED_TIME;
 extern bool GNSS_DENIED_FLAG;
 extern double SYS_STOP_TIME;
 extern double VINS_INIT_TIME;
 extern double VINS_INIT_TIMESTAMP;
 extern bool READ_CSV_FILE;
 extern std::string CSV_FILENAME;
 
 void readParameters(ros::NodeHandle &n);
 
 enum SIZE_PARAMETERIZATION
 {
     SIZE_POSE = 7,
     SIZE_SPEEDBIAS = 9,
     SIZE_FEATURE = 1
 };
 
 enum StateOrder
 {
     O_P = 0,
     O_R = 3,
     O_V = 6,
     O_BA = 9,
     O_BG = 12
 };
 
 enum NoiseOrder
 {
     O_AN = 0,
     O_GN = 3,
     O_AW = 6,
     O_GW = 9
 };