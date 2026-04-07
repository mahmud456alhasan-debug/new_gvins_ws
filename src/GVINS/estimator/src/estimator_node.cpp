#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include <gnss_comm/gnss_ros.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <gvins/LocalSensorExternalTrigger.h>
#include <sensor_msgs/NavSatFix.h>
#include <gnss_comm/GnssPVTSolnMsg.h>

#include "estimator.h"
#include "parameters.h"
#include "utility/visualization.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>  // for std::setprecision
#include <cstdlib>
#include <chrono>

// 定义一个结构体来存储每行数据
struct DataRow {
    int frame_count;
    double gnss_time;
    double px, py, pz;
    double vx, vy, vz;
    double yaw, pitch, roll;
    double bax, bay, baz;
    double bgx, bgy, bgz;
};

// CSV realted
std::ifstream csv_file;

using namespace gnss_comm;

#define MAX_GNSS_CAMERA_DELAY 0.05 // 0.05~20Hz  0.033~30Hz 

std::unique_ptr<Estimator> estimator_ptr;

std::condition_variable con;
double current_time = -1;
queue<sensor_msgs::ImuConstPtr> imu_buf;
queue<sensor_msgs::PointCloudConstPtr> feature_buf;
queue<std::vector<ObsPtr>> gnss_meas_buf;
queue<sensor_msgs::PointCloudConstPtr> relo_buf;
int sum_of_wait = 0;

std::mutex m_buf;
std::mutex m_state;
std::mutex i_buf;
std::mutex m_estimator;

double latest_time;
Eigen::Vector3d tmp_P;
Eigen::Quaterniond tmp_Q;
Eigen::Vector3d tmp_V;
Eigen::Vector3d tmp_Ba;
Eigen::Vector3d tmp_Bg;
Eigen::Vector3d acc_0;
Eigen::Vector3d gyr_0;
bool init_feature = 0;
bool init_imu = 1;
double last_imu_t = -1;

std::mutex m_time;
double next_pulse_time;
bool next_pulse_time_valid;
double time_diff_gnss_local;
bool time_diff_valid;
double latest_gnss_time;
double tmp_last_feature_time;
uint64_t feature_msg_counter;
int skip_parameter;

static void pvt_adaptive_callback(const gnss_comm::GnssPVTSolnMsgConstPtr &msg)
{
    if (!GNSS_ADAPTIVE_CARR_SOLN_SWITCH)
        return;
    const bool want_tight = (msg->carr_soln < 2);
    const bool prev = GNSS_ENABLE.exchange(want_tight);
    if (prev != want_tight)
        ROS_WARN("GNSS adaptive: carr_soln=%u -> raw GNSS in estimator %s", (unsigned)msg->carr_soln, want_tight ? "ON" : "OFF");
}

void write_imu_pre_state_ai(const long long int& total_frame_cout, const double& time, const Vector3d& pos, const Vector3d& vel, const Matrix3d& rot_mat, const double& yaw_enu_local) {
    static bool is_has_header = false;
    ofstream of_file(RECORD_FILE_PATH + "/imu_pre_state_ai.csv", ios::app);
    
    Eigen::Vector3d temp_eulerAngle = rot_mat.eulerAngles(2, 1, 0); // ZYX 顺序，即 yaw pitch roll
    
    if(of_file.is_open()) {
        if(!is_has_header) {
            of_file  << "frame_count" << ","
                    <<  "gnss_time"<< ","
                    <<  "px"  << ","
                    <<  "py"  << ","
                    <<  "pz"  << ","
                    <<  "yaw_local2enu"  << ","
                    <<  "vx"  << ","
                    <<  "vy"  << ","
                    <<  "vz"  << ","
                    <<  "pitch_local"  << ","
                    <<  "roll_local"  << ","
                    <<  "yaw_local"  << std::endl;
            is_has_header = true;
        } else {
            of_file  << std::setprecision(16)   
                        << total_frame_cout << ","
                        << time << ","
                        << pos.x() << ","
                        << pos.y() << ","
                        << pos.z() << ","
                        << yaw_enu_local << ","
                        << vel.x() << ","
                        << vel.y() << ","
                        << vel.z() << ","
                        << temp_eulerAngle.x() << ","
                        << temp_eulerAngle.y() << ","
                        << temp_eulerAngle.z() << std::endl;
                        
        }
        of_file.close();
    } else {
            std::cerr << "Failed to open imu_pre_state_ai.csv!" << std::endl;
    }
}

// 保存 PVQ 状态
void write_PVQ_state(string filename, const long long int& total_frame_cout, const Estimator::SolverFlag& solver_flag, const double& time,
                        const Vector3d& pos, const Matrix3d& rot, const Vector3d& vel) {
    static bool is_has_header = false;
    
    std::ofstream of_imu_preitgr(RECORD_FILE_PATH + filename, std::ios::app);
    
    Eigen::Quaterniond rot_quater(rot);
    
    if(of_imu_preitgr.is_open()) {
        if(!is_has_header) {
            of_imu_preitgr  << "frame_count" << ","
                    << "sys_state"   << ","
                    << "time" << ","
                    <<  "px"  << ","
                    <<  "py"  << ","
                    <<  "pz"  << ","
                    <<  "qx"  << ","
                    <<  "qy"  << ","
                    <<  "qz"  << ","
                    <<  "qw"  << ","
                    <<  "vx"  << ","
                    <<  "vy"  << ","
                    <<  "vz"  << std::endl;
            is_has_header = true;
        } else {
            of_imu_preitgr  << std::setprecision(16)
                    << total_frame_cout << ","
                    << solver_flag << ","
                    << time  << ","
                    << pos.x()  << ","
                    << pos.y()  << ","
                    << pos.z()  << ","
                    << rot_quater.coeffs().x() << ","
                    << rot_quater.coeffs().y() << ","
                    << rot_quater.coeffs().z() << ","
                    << rot_quater.coeffs().w() << ","
                    << vel.x()  << ","
                    << vel.y()  << ","
                    << vel.z()  << std::endl;
            }
            of_imu_preitgr.close();
        } else {
            std::cerr << "Failed to open" << filename << std::endl;
    }
    
}

void predict(const sensor_msgs::ImuConstPtr &imu_msg)
{
    double t = imu_msg->header.stamp.toSec();
    if (init_imu)
    {
        latest_time = t;
        init_imu = 0;
        return;
    }
    double dt = t - latest_time;
    latest_time = t;

    double dx = imu_msg->linear_acceleration.x;
    double dy = imu_msg->linear_acceleration.y;
    double dz = imu_msg->linear_acceleration.z;
    Eigen::Vector3d linear_acceleration{dx, dy, dz};

    double rx = imu_msg->angular_velocity.x;
    double ry = imu_msg->angular_velocity.y;
    double rz = imu_msg->angular_velocity.z;
    Eigen::Vector3d angular_velocity{rx, ry, rz};

    Eigen::Vector3d un_acc_0 = tmp_Q * (acc_0 - tmp_Ba) - estimator_ptr->g;

    Eigen::Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - tmp_Bg;
    tmp_Q = tmp_Q * Utility::deltaQ(un_gyr * dt);

    Eigen::Vector3d un_acc_1 = tmp_Q * (linear_acceleration - tmp_Ba) - estimator_ptr->g;

    Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);

    tmp_P = tmp_P + dt * tmp_V + 0.5 * dt * dt * un_acc;
    tmp_V = tmp_V + dt * un_acc;

    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
}

void update()
{
    TicToc t_predict;
    latest_time = current_time;
    tmp_P = estimator_ptr->Ps[WINDOW_SIZE];
    tmp_Q = estimator_ptr->Rs[WINDOW_SIZE];
    tmp_V = estimator_ptr->Vs[WINDOW_SIZE];
    tmp_Ba = estimator_ptr->Bas[WINDOW_SIZE];
    tmp_Bg = estimator_ptr->Bgs[WINDOW_SIZE];
    acc_0 = estimator_ptr->acc_0;
    gyr_0 = estimator_ptr->gyr_0;

    queue<sensor_msgs::ImuConstPtr> tmp_imu_buf = imu_buf;
    for (sensor_msgs::ImuConstPtr tmp_imu_msg; !tmp_imu_buf.empty(); tmp_imu_buf.pop())
        predict(tmp_imu_buf.front());

}

bool getMeasurements(std::vector<sensor_msgs::ImuConstPtr> &imu_msg, sensor_msgs::PointCloudConstPtr &img_msg, std::vector<ObsPtr> &gnss_msg)
{
    // 1. 队列空判定
    if (imu_buf.empty() || feature_buf.empty() || (GNSS_ENABLE.load() && gnss_meas_buf.empty())) // 当IMU、图像和GNSS中只要有一个数据缓存器为空，那么直接返回false
        return false;
    
    // 2. 对齐IMU和image时间
    double front_feature_ts = feature_buf.front()->header.stamp.toSec(); // feature_buf 这个 queue 中第一帧特征点的时间戳

    if (!(imu_buf.back()->header.stamp.toSec() > front_feature_ts)) // 确保queue中最后一帧IMU数据的时间要晚于 第一帧特征点的时间戳，否则就等待IMU
    {
        //ROS_WARN("wait for imu, only should happen at the beginning");
        sum_of_wait++;
        return false;
    }
    double front_imu_ts = imu_buf.front()->header.stamp.toSec(); // queue中第一帧imu时间
    while (!feature_buf.empty() && front_imu_ts > front_feature_ts) // 图像queue不为空 但是 第一帧IMU时间 大于 图像的第一帧时间，此时需要丢弃部分图像
    {
        ROS_WARN("throw img, only should happen at the beginning");
        feature_buf.pop();
        front_feature_ts = feature_buf.front()->header.stamp.toSec();
    }

    // 3. 对齐GNSS和image , 进而实现 GNSS、IMU、image三者对齐
    if (GNSS_ENABLE.load())
    {
        front_feature_ts += time_diff_gnss_local; // 修正图像帧的时间戳 （另一种说法是：UTC和GPS之间相差了18.0s，需要看看ROS读取的时间是不是和GPS时间相差了18.0s，有可能不差）
        double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); 
        while (!gnss_meas_buf.empty() && front_gnss_ts < front_feature_ts-MAX_GNSS_CAMERA_DELAY) //直接考虑GNSS和相机触发时间之间最大的延时MAX_GNSS_CAMERA_DELAY
        { // 如果front_gnss_ts < front_feature_ts-MAX_GNSS_CAMERA_DELAY，那么说明GNSS中包含过旧的数据，直接丢弃
            // ROS_WARN("throw gnss, only should happen at the beginning");
            gnss_meas_buf.pop();
            if (gnss_meas_buf.empty()) return false;
            front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time);
        }
        if (gnss_meas_buf.empty()) // gnss数据空判断，这行代码可能没啥用？
        {
            ROS_WARN("wait for gnss...");
            return false;
        }
        else if (abs(front_gnss_ts-front_feature_ts) < MAX_GNSS_CAMERA_DELAY)
        { // 如果GNSS第一个数据的时间戳和修正后的第一帧图像时间戳之间的差之在MAX_GNSS_CAMERA_DELAY范围内，不满足上述while的循环条件，所以无法丢弃
          // 保存第一个GNSS数据，并将其从queue中pop出
            gnss_msg = gnss_meas_buf.front(); //提取GNSS数据
            gnss_meas_buf.pop();
        }
    }

    img_msg = feature_buf.front(); //提取 feature数据
    feature_buf.pop();

    // 提取两帧image之间的imu数据，便于进行预积分
    while (imu_buf.front()->header.stamp.toSec() < img_msg->header.stamp.toSec() + estimator_ptr->td)
    {
        imu_msg.emplace_back(imu_buf.front());
        imu_buf.pop();
    }
    imu_msg.emplace_back(imu_buf.front());
    if (imu_msg.empty())
        ROS_WARN("no imu between two image");
    return true;
}

void imu_callback(const sensor_msgs::ImuConstPtr &imu_msg)
{
    if (imu_msg->header.stamp.toSec() <= last_imu_t)
    {
        ROS_WARN("imu message in disorder!");
        return;
    }
    last_imu_t = imu_msg->header.stamp.toSec();

    m_buf.lock();
    imu_buf.push(imu_msg);
    m_buf.unlock();
    con.notify_one();

    last_imu_t = imu_msg->header.stamp.toSec();

    {
        std::lock_guard<std::mutex> lg(m_state);
        predict(imu_msg);
        std_msgs::Header header = imu_msg->header;
        header.frame_id = "world";
        if (estimator_ptr->solver_flag == Estimator::SolverFlag::NON_LINEAR)
            pubLatestOdometry(tmp_P, tmp_Q, tmp_V, header);
    }
}

// GPS BEIDOU 星历回调
void gnss_ephem_callback(const GnssEphemMsgConstPtr &ephem_msg)
{
    EphemPtr ephem = msg2ephem(ephem_msg);
    estimator_ptr->inputEphem(ephem);
}

// 格罗纳斯 星历回调
void gnss_glo_ephem_callback(const GnssGloEphemMsgConstPtr &glo_ephem_msg)
{
    GloEphemPtr glo_ephem = msg2glo_ephem(glo_ephem_msg);
    estimator_ptr->inputEphem(glo_ephem);
}

// 电离层参数回调
void gnss_iono_params_callback(const StampedFloat64ArrayConstPtr &iono_msg)
{
    double ts = iono_msg->header.stamp.toSec();
    std::vector<double> iono_params;
    std::copy(iono_msg->data.begin(), iono_msg->data.end(), std::back_inserter(iono_params));
    assert(iono_params.size() == 8);
    estimator_ptr->inputIonoParams(ts, iono_params);
}

// 原始量测回调
void gnss_meas_callback(const GnssMeasMsgConstPtr &meas_msg)
{
    std::vector<ObsPtr> gnss_meas = msg2meas(meas_msg); // 从ROS消息解析GNSS测量值

    latest_gnss_time = time2sec(gnss_meas[0]->time);

    // cerr << "gnss ts is " << std::setprecision(20) << time2sec(gnss_meas[0]->time) << endl;
    if (!time_diff_valid)   return;

    m_buf.lock();
    gnss_meas_buf.push(std::move(gnss_meas)); //
    m_buf.unlock();
    con.notify_one(); // 唤醒一个正在等待的线程，使其获得锁并继续执行
}

void feature_callback(const sensor_msgs::PointCloudConstPtr &feature_msg)
{
    ++ feature_msg_counter;

    if (skip_parameter < 0 && time_diff_valid)
    {
        const double this_feature_ts = feature_msg->header.stamp.toSec()+time_diff_gnss_local;
        if (latest_gnss_time > 0 && tmp_last_feature_time > 0)
        {
            if (abs(this_feature_ts - latest_gnss_time) > abs(tmp_last_feature_time - latest_gnss_time))
                skip_parameter = feature_msg_counter%2;       // skip this frame and afterwards
            else
                skip_parameter = 1 - (feature_msg_counter%2);   // skip next frame and afterwards
        }
        // cerr << "feature counter is " << feature_msg_counter << ", skip parameter is " << int(skip_parameter) << endl;
        tmp_last_feature_time = this_feature_ts;
    }

    if (skip_parameter >= 0 && int(feature_msg_counter%2) != skip_parameter)
    {
        m_buf.lock();
        feature_buf.push(feature_msg);
        m_buf.unlock();
        con.notify_one();
    }
}

void local_trigger_info_callback(const gvins::LocalSensorExternalTriggerConstPtr &trigger_msg)
{// trigger_msg记录了相机被GNSS的脉冲触发的时间，也可以理解成图像的命名，由于存在硬件延时，和真正的GNSS时间有区别
 // 这也是为什么后面校正local和world时间的原因
    std::lock_guard<std::mutex> lg(m_time);

    if (next_pulse_time_valid)
    {
        // 计算local和gnss的时间差
        time_diff_gnss_local = next_pulse_time - trigger_msg->header.stamp.toSec();
        estimator_ptr->inputGNSSTimeDiff(time_diff_gnss_local);
        if (!time_diff_valid)       // just get calibrated
            std::cout << "time difference between GNSS and VI-Sensor got calibrated: "
                << std::setprecision(15) << time_diff_gnss_local << " s\n";
        time_diff_valid = true;
    }
}

void gnss_tp_info_callback(const GnssTimePulseInfoMsgConstPtr &tp_msg)
{
    gtime_t tp_time = gpst2time(tp_msg->time.week, tp_msg->time.tow);
    if (tp_msg->utc_based || tp_msg->time_sys == SYS_GLO)
        tp_time = utc2gpst(tp_time);
    else if (tp_msg->time_sys == SYS_GAL)
        tp_time = gst2time(tp_msg->time.week, tp_msg->time.tow);
    else if (tp_msg->time_sys == SYS_BDS)
        tp_time = bdt2time(tp_msg->time.week, tp_msg->time.tow);
    else if (tp_msg->time_sys == SYS_NONE)
    {
        std::cerr << "Unknown time system in GNSSTimePulseInfoMsg.\n";
        return;
    }
    double gnss_ts = time2sec(tp_time);

    std::lock_guard<std::mutex> lg(m_time);
    next_pulse_time = gnss_ts;
    next_pulse_time_valid = true;
}

void restart_callback(const std_msgs::BoolConstPtr &restart_msg)
{
    if (restart_msg->data == true)
    {
        ROS_WARN("restart the estimator!");
        m_buf.lock();
        while(!feature_buf.empty())
            feature_buf.pop();
        while(!imu_buf.empty())
            imu_buf.pop();
        m_buf.unlock();
        m_estimator.lock();
        estimator_ptr->clearState();
        estimator_ptr->setParameter();
        m_estimator.unlock();
        current_time = -1;
        last_imu_t = 0;
    }
    return;
}


void process()
{
    while (true) // 这个循环相当于 measurement_process 线程的while(1)
    {
        std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr>> measurements; //用于将一[组]IMU消息和一条点云(或图像？)消息关联
        std::vector<sensor_msgs::ImuConstPtr> imu_msg; //用于存储一组[IMU]数据
        sensor_msgs::PointCloudConstPtr img_msg; 
        std::vector<ObsPtr> gnss_msg; // GNSS观测数据
                
        // 1. IMU GNSS image 数据准备
        // 获取三者观测数据，为融合做准备，只要数据不满足要求，就一直在这里等待
        // 这段代码的作用是在互斥锁的保护下等待，直到满足获取测量数据的条件为止，然后解锁互斥锁以允许其他线程访问测量数据。这是多线程编程中常见的同步机制，以确保数据的安全访问。
        std::unique_lock<std::mutex> lk(m_buf); // unique_lock是一个在互斥锁上的RAII（资源获取即初始化）封装，用于在作用域结束时自动释放锁
        con.wait(lk, [&]
                 {
                    return getMeasurements(imu_msg, img_msg, gnss_msg); // 根据时间戳判断量测信息的合法性
                 }); // 使用条件变量con等待条件成立， 它在lk锁的保护下等待，直到提供的lambda函数返回true; 这个等待的目的是为了等待测量数据就绪，以便后续代码可以安全的访问这些数据
        lk.unlock(); //解锁lk

        m_estimator.lock();
        
        // 2. IMU数据处理 current_time是目前算法处理到的IMU的位置
        // 根据IMU数据的时间戳，判断是进行IMU数据的处理还是进行线性插值
        double dx = 0, dy = 0, dz = 0, rx = 0, ry = 0, rz = 0; // 初始化线性加速度和角速度的变量。
        for (auto &imu_data : imu_msg)
        {
            double t = imu_data->header.stamp.toSec(); // 获取当前IMU消息的时间戳。
            double img_t = img_msg->header.stamp.toSec() + estimator_ptr->td; //计算图像消息的时间戳，考虑了一个时间偏移td
            if (t <= img_t)  // 若IMU时间滞后于image时间，则直接对IMU进行预积分
            { 
                if (current_time < 0) // 如果当前时间current_time小于0，说明是第一次处理IMU数据，将当前时间设置为当前IMU消息的时间戳。
                    current_time = t;
                double dt = t - current_time; // 计算两次IMU消息之间的时间差
                ROS_ASSERT(dt >= 0); // 断言确保时间差为非负数，即时间戳是递增的
                current_time = t; // 更新当前IMU预积分处理到的位置
                dx = imu_data->linear_acceleration.x; // 依次取出数据并进行变量赋值
                dy = imu_data->linear_acceleration.y;
                dz = imu_data->linear_acceleration.z;
                rx = imu_data->angular_velocity.x;
                ry = imu_data->angular_velocity.y;
                rz = imu_data->angular_velocity.z;
                estimator_ptr->processIMU(dt, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz));
                //printf("imu: dt:%f a: %f %f %f w: %f %f %f\n",dt, dx, dy, dz, rx, ry, rz);

            }
            else // 如果IMU时间超前于image时间，
            {
                double dt_1 = img_t - current_time; // 预积分的前半段 dt1：本次图像时间戳 - 上次IMU预积分到的位置（上次预积分的终点保存在current时间），这部分是从[上次IMU预积分终点]到[本次image图像时间戳]的位置
                double dt_2 = t - img_t; // 预积分的后半段 dt2 : 本次IMU时间戳 - 本次图像时间戳
                current_time = img_t; 
                ROS_ASSERT(dt_1 >= 0);
                ROS_ASSERT(dt_2 >= 0);
                ROS_ASSERT(dt_1 + dt_2 > 0);
                double w1 = dt_2 / (dt_1 + dt_2);
                double w2 = dt_1 / (dt_1 + dt_2);
                // 对线性加速度进行插值。新的线性加速度是前一时刻线性加速度与当前IMU消息线性加速度的线性组合。
                // TODO 没理解这部分为什么这么处理
                dx = w1 * dx + w2 * imu_data->linear_acceleration.x;
                dy = w1 * dy + w2 * imu_data->linear_acceleration.y;
                dz = w1 * dz + w2 * imu_data->linear_acceleration.z;
                rx = w1 * rx + w2 * imu_data->angular_velocity.x;
                ry = w1 * ry + w2 * imu_data->angular_velocity.y;
                rz = w1 * rz + w2 * imu_data->angular_velocity.z;
                
                estimator_ptr->processIMU(dt_1, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz)); //IMU预积分
                //printf("dimu: dt:%f a: %f %f %f w: %f %f %f\n",dt_1, dx, dy, dz, rx, ry, rz);
            }
        }

        // 3. 处理GNSS数据，使得GNSS数据可用 
        if (GNSS_ENABLE.load() && !gnss_msg.empty())
            estimator_ptr->processGNSS(gnss_msg);

        // 4. 处理image数据
        // 遍历图像点云数据，提取特征点信息，组织成 image 数据结构，然后调用 estimator_ptr->processImage(image, img_msg->header); 进行状态估计
        ROS_DEBUG("processing vision data with stamp %f \n", img_msg->header.stamp.toSec());

        TicToc t_s; // 计时器类，通过记录开始和结束的时间点来测量代码段的执行时间。
        map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> image; // 使用一个map数据结构，其中键是特征ID，值是一个向量，
                                                                        // 包含相机ID和一个7x1的矩阵，表示特征的空间坐标（x, y, z），图像坐标（p_u, p_v），以及速度信息（velocity_x, velocity_y）。
        
        for (unsigned int i = 0; i < img_msg->points.size(); i++) // 遍历图像消息中的所有点。
        {
            int v = img_msg->channels[0].values[i] + 0.5; // 从通道中获取特征ID（假设通道0中存储了特征ID），并进行四舍五入操作。
            int feature_id = v / NUM_OF_CAM; // 从特征ID中提取特征在相机中的ID
            int camera_id = v % NUM_OF_CAM; // 相机ID
            double x = img_msg->points[i].x; // 获取特征点的空间坐标。
            double y = img_msg->points[i].y;
            double z = img_msg->points[i].z;
            double p_u = img_msg->channels[1].values[i]; // 获取特征点在图像上的坐标。
            double p_v = img_msg->channels[2].values[i];
            double velocity_x = img_msg->channels[3].values[i]; // 获取特征点的速度信息。
            double velocity_y = img_msg->channels[4].values[i];
            ROS_ASSERT(z == 1); // 进行断言检查，确保特征点的深度（z坐标）为1。
            Eigen::Matrix<double, 7, 1> xyz_uv_velocity; // 将特征点的信息组合成一个7x1的矩阵。
            xyz_uv_velocity << x, y, z, p_u, p_v, velocity_x, velocity_y; 
            image[feature_id].emplace_back(camera_id,  xyz_uv_velocity); // 将特征点的信息存储到 image 中，以特征ID为键，值为包含相机ID和特征信息的向量。
        }
        estimator_ptr->processImage(image, img_msg->header); // 将整理好的特征数据传递给估计器的 processImage 函数进行处理。
                                                             // 这个函数用于将视觉信息与IMU数据进行融合，用于进行视觉-惯性里程计（Visual-Inertial Odometry，VIO）或者类似的定位与导航任务。

        // 5. 发布消息，用于rviz显示
        double whole_t = t_s.toc();
        printStatistics(*estimator_ptr, whole_t); // 数据静态打印（调试用）
        std_msgs::Header header = img_msg->header;
        header.frame_id = "world";
        
        pubOdometry(*estimator_ptr, header);
        pubGnssResult(*estimator_ptr, header);
        pubKeyPoses(*estimator_ptr, header);
        pubCameraPose(*estimator_ptr, header);
        pubPointCloud(*estimator_ptr, header);
        pubTF(*estimator_ptr, header);
        pubKeyframe(*estimator_ptr);

        estimator_ptr->total_frame_count = estimator_ptr->total_frame_count + 1;

        std::vector<DataRow> temp;

        m_estimator.unlock();
        m_buf.lock();
        m_state.lock();

        // 6. 如果估计器的求解标志为非线性，调用 update() 函数进行状态更新
        if (estimator_ptr->solver_flag == Estimator::SolverFlag::NON_LINEAR) 
            update();
        m_state.unlock();
        m_buf.unlock();
    }
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "gvins");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);
    readParameters(n); // 从ROS参数服务器中读物参数，这里的大部分参数都来自 config文件夹下的yaml文件
    estimator_ptr.reset(new Estimator());
    estimator_ptr->setParameter();
#ifdef EIGEN_DONT_PARALLELIZE // eigen 并行化运算判断
    ROS_DEBUG("EIGEN_DONT_PARALLELIZE");
#endif

    registerPub(n); // 注册要发布的话题，详细的定义在visualization.cpp中

    next_pulse_time_valid = false; //GNSS时间脉冲是否有效的标志位
    time_diff_valid = false;  //GNSS时间与本地（固定）时间差（原文是18.0s）是否有效
    latest_gnss_time = -1; //最新一帧的GNSS时间
    tmp_last_feature_time = -1; //TODO: tmp_last_feature_time 啥意思
    feature_msg_counter = 0; //前端特征点帧数计数，这意味着前端每来一帧新的特征点，该计数值+1

    if (GNSS_ENABLE.load() || GNSS_ADAPTIVE_CARR_SOLN_SWITCH)
        skip_parameter = -1; // TODO：skip_parameter啥意思
    else
        skip_parameter = 0;


    // subsrciber参数详解：topic-订阅的节点名；queue_size-待处理信息队列大小；callback-回调函数
    ros::Subscriber sub_imu = n.subscribe(IMU_TOPIC, 2000, imu_callback, ros::TransportHints().tcpNoDelay()); //订阅IMU信息,TransportHints().tcpNoDelay() 是用于设置消息传输选项的一种方式，这里是设置消息传输时不使用TCP延迟。
    ros::Subscriber sub_feature = n.subscribe("/gvins_feature_tracker/feature", 2000, feature_callback); // 订阅前端特征点
    ros::Subscriber sub_restart = n.subscribe("/gvins_feature_tracker/restart", 2000, restart_callback); // 订阅故障信息

    // 订阅GNSS信息
    ros::Subscriber sub_ephem, sub_glo_ephem, sub_gnss_meas, sub_gnss_iono_params;
    ros::Subscriber sub_gnss_time_pluse_info, sub_local_trigger_info;
    ros::Subscriber sub_pvt;
    if (GNSS_ADAPTIVE_CARR_SOLN_SWITCH)
        sub_pvt = n.subscribe(GNSS_PVT_TOPIC, 500, pvt_adaptive_callback, ros::TransportHints().tcpNoDelay());

    if (GNSS_ENABLE.load() || GNSS_ADAPTIVE_CARR_SOLN_SWITCH)
    {
        // 订阅两个不同的星历话题，是因为两个导航系统下的星历格式不一样
        // 订阅星历信息：卫星的位置、速度、时间偏差等信息
        sub_ephem = n.subscribe(GNSS_EPHEM_TOPIC, 100, gnss_ephem_callback);  // GNSS星历信息
        sub_glo_ephem = n.subscribe(GNSS_GLO_EPHEM_TOPIC, 100, gnss_glo_ephem_callback); // GLO：GLONASS。格洛纳斯星历信息
        sub_gnss_meas = n.subscribe(GNSS_MEAS_TOPIC, 100, gnss_meas_callback); // 卫星的观测信息 raw measurement
        sub_gnss_iono_params = n.subscribe(GNSS_IONO_PARAMS_TOPIC, 100, gnss_iono_params_callback);  // 电离层参数订阅

        // GNSS和VIO的时间是否同步判断
        if (GNSS_LOCAL_ONLINE_SYNC) // 在线同步
        {
            sub_gnss_time_pluse_info = n.subscribe(GNSS_TP_INFO_TOPIC, 100, gnss_tp_info_callback); // 订阅GNSS脉冲信息
            sub_local_trigger_info = n.subscribe(LOCAL_TRIGGER_INFO_TOPIC, 100, local_trigger_info_callback); // 订阅相机触发时间
        }
        else
        {
            time_diff_gnss_local = GNSS_LOCAL_TIME_DIFF; // 18.0s的差值
            estimator_ptr->inputGNSSTimeDiff(time_diff_gnss_local); // 本函数直接将 diff_t_gnss_local 写死为 GNSS_LOCAL_TIME_DIFF（18.0s）
            time_diff_valid = true;
        }
        // annotation: GVINS 运行时使用了GNSS和VIO的结果，但是两者是不同空间的产物，必然会存在时间差，需要进行补偿

    }


    std::thread measurement_process{process}; // measurement_process是被创建的线程名，process是要被执行的函数

    ros::spin(); //阻塞主线程，等待处理回调函数

    return 0;
}