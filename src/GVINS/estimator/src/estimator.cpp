#include "estimator.h"

Estimator::Estimator(): f_manager{Rs}
{
    ROS_INFO("init begins");
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
        pre_integrations[i] = nullptr;
    clearState(); //将所有的状态置0
    ROS_INFO("[Damon] All states are set to zero!");
}

void Estimator::setParameter()
{
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = TIC[i];
        ric[i] = RIC[i];
    }
    f_manager.setRic(ric);
    ProjectionFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    ProjectionTdFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    td = TD;
}

void Estimator::clearState()
{
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        Rs[i].setIdentity();
        Ps[i].setZero();
        Vs[i].setZero();
        Bas[i].setZero();
        Bgs[i].setZero();
        dt_buf[i].clear();
        linear_acceleration_buf[i].clear();
        angular_velocity_buf[i].clear();

        if (pre_integrations[i] != nullptr)
            delete pre_integrations[i];
        pre_integrations[i] = nullptr;
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d::Zero();
        ric[i] = Matrix3d::Identity();
    }

    for (auto &it : all_image_frame)
    {
        if (it.second.pre_integration != nullptr)
        {
            delete it.second.pre_integration;
            it.second.pre_integration = nullptr;
        }
    }

    solver_flag = INITIAL;
    first_imu = false,
    sum_of_back = 0;
    sum_of_front = 0;
    frame_count = 0;
    solver_flag = INITIAL;
    initial_timestamp = 0;
    all_image_frame.clear();
    td = TD;

    gnss_ready = false;
    anc_ecef.setZero();
    R_ecef_enu.setIdentity();
    para_yaw_enu_local[0] = 0;
    yaw_enu_local = 0;
    sat2ephem.clear();
    sat2time_index.clear();
    sat_track_status.clear();
    latest_gnss_iono_params.clear();
    std::copy(GNSS_IONO_DEFAULT_PARAMS.begin(), GNSS_IONO_DEFAULT_PARAMS.end(), 
        std::back_inserter(latest_gnss_iono_params));
    diff_t_gnss_local = 0;

    first_optimization = true;

    if (tmp_pre_integration != nullptr)
        delete tmp_pre_integration;
    if (last_marginalization_info != nullptr)
        delete last_marginalization_info;

    tmp_pre_integration = nullptr;
    last_marginalization_info = nullptr;
    last_marginalization_parameter_blocks.clear();

    f_manager.clearState();

    failure_occur = 0;
}

// IMU预积分
void Estimator::processIMU(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{
    if (!first_imu)  // 检查是否是第一次接收到IMU数据。如果是第一次，进行初始化，记录当前的线性加速度 acc_0 和角速度 gyr_0
    {
        first_imu = true;
        acc_0 = linear_acceleration;
        gyr_0 = angular_velocity;
    }

    if (!pre_integrations[frame_count]) // 检查是否已经为当前帧创建了 pre_integrations 对象，用于存储积分信息。如果没有，则创建一个新的 IntegrationBase 对象。
    {
        pre_integrations[frame_count] = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};
    }
    if (frame_count != 0) // 除第一帧之外的情况（因为第一帧没有前一帧的状态作为参考）
    {
        pre_integrations[frame_count]->push_back(dt, linear_acceleration, angular_velocity); // 将加速度、角速度进行预积分。注意：此处的push_back非C++ lib中的push_back，
                                                                                             // 这里pushback之后已经进行了预积分
        //if(solver_flag != NON_LINEAR)
            tmp_pre_integration->push_back(dt, linear_acceleration, angular_velocity); // 临时的预积分器，这个积分器是给processImage函数使用的

        // 将当前帧的【时间间隔】、【线性加速度】和【角速度】添加到相应的缓冲区中，用于后续的处理。
        dt_buf[frame_count].push_back(dt);
        linear_acceleration_buf[frame_count].push_back(linear_acceleration);
        angular_velocity_buf[frame_count].push_back(angular_velocity);

        int j = frame_count;         
        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - g; // 将前一帧的线性加速度 acc_0 通过旋转矩阵 Rs[j] 变换到当前帧坐标系，并减去当前帧的加速度偏置 Bas[j] 和重力向量 g。
        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j]; // 将前一帧的角速度 gyr_0 和当前帧的角速度取平均，并减去当前帧的陀螺仪偏置 Bgs[j]。
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix(); // 更新当前帧的旋转矩阵 Rs[j]，这里使用了一个旋转增量 Utility::deltaQ(un_gyr * dt)。
        Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - g; // 将当前帧的线性加速度通过旋转矩阵 Rs[j] 变换到当前帧坐标系，并减去当前帧的加速度偏置 Bas[j] 和重力向量 g。
        
        // 状态量更新
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1); // 加速度更新：取前一帧和当前帧的变换后的加速度的平均值，得到当前帧的变换后的加速度
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc; // 位置更新： 根据当前帧的速度 Vs[j] 和变换后的加速度 un_acc 更新当前帧的位置 Ps[j]
        Vs[j] += dt * un_acc; // 速度更新：根据变换后的加速度 un_acc 更新当前帧的速度 Vs[j]
    }
    // 更新前一帧的线性加速度和角速度为当前帧的值，以备下一次迭代使用。
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
}

// 实现视觉惯性里程计（Visual-Inertial Odometry，VIO）中的图像处理和状态估计的一系列步骤
void Estimator::processImage(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, const std_msgs::Header &header)
{
    ROS_DEBUG("new image coming ------------------------------------------");
    ROS_DEBUG("Adding feature points %lu", image.size());
    if (f_manager.addFeatureCheckParallax(frame_count, image, td)) // 检查新图像中的特征点，并检查其视差是否足够，以决定是否进行边缘化（marginalization）。
        marginalization_flag = MARGIN_OLD; // 如果返回值为 true，说明成功添加特征并通过视差检查，当前帧满足条件，可以边缘化掉一些旧的状态，减少计算开销
    else
        marginalization_flag = MARGIN_SECOND_NEW; // 如果返回值为 false，说明未通过视差检查，当前帧的视差不够，不能边缘化旧状态，需要保留更多的新状态进行后续计算

    ROS_DEBUG("this frame is--------------------%s", marginalization_flag ? "reject" : "accept");
    ROS_DEBUG("%s", marginalization_flag ? "Non-keyframe" : "Keyframe");
    ROS_DEBUG("Solving %d", frame_count);
    ROS_DEBUG("number of feature: %d", f_manager.getFeatureCount());
    Headers[frame_count] = header; // 将当前帧的头信息保存到全局变量中

    ImageFrame imageframe(image, header.stamp.toSec()); // 创建一个 ImageFrame 对象，其中包含当前图像的特征信息和时间戳。
    imageframe.pre_integration = tmp_pre_integration; // 将临时预积分数据（tmp_pre_integration）赋值给 ImageFrame 对象。
    all_image_frame.insert(make_pair(header.stamp.toSec(), imageframe)); // 将 ImageFrame 对象插入到一个按时间戳排序的容器中，以备后续使用。
    tmp_pre_integration = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]}; // 为下一帧创建新的预积分对象。

    if(ESTIMATE_EXTRINSIC == 2) // 检查是否需要校准外部参数，其中 ESTIMATE_EXTRINSIC 的值为2表示需要进行外部参数的校准。
    {
        ROS_INFO("calibrating extrinsic param, rotation movement is needed");
        if (frame_count != 0) // 确保不是第一帧。
        {
            vector<pair<Vector3d, Vector3d>> corres = f_manager.getCorresponding(frame_count - 1, frame_count); // 获取前一帧和当前帧之间的对应特征点。
            Matrix3d calib_ric; // 定义用于校准的旋转矩阵
            if (initial_ex_rotation.CalibrationExRotation(corres, pre_integrations[frame_count]->delta_q, calib_ric)) // 调用外部旋转校准函数。
            {
                // 如果校准成功，更新外部旋转矩阵，并将 ESTIMATE_EXTRINSIC 设置为1。
                ROS_WARN("initial extrinsic rotation calib success");
                ROS_WARN_STREAM("initial extrinsic rotation: " << endl << calib_ric);
                ric[0] = calib_ric;
                RIC[0] = calib_ric;
                ESTIMATE_EXTRINSIC = 1;
            }
        }
    }

    if (solver_flag == INITIAL) // 检查是否处于初始化阶段。
    {
        if (frame_count == WINDOW_SIZE) // 检查是否已经收集了足够的帧用于初始化。
        {
            bool result = false; // 初始化结果标志。
            if( ESTIMATE_EXTRINSIC != 2 && (header.stamp.toSec() - initial_timestamp) > 0.1) // 检查是否需要进行外部参数校准，以及是否已经过了一定时间。
            {
                result = initialStructure(); // 调用初始化结构的函数。
                initial_timestamp = header.stamp.toSec(); // 更新初始化时间戳
            }
            
            if(result) // 如果初始化成功
            {
                solver_flag = NON_LINEAR;  //设置求解标志为非线性求解。 
                solveOdometry(); 
                slideWindow();
                f_manager.removeFailures(); // 移除特征点管理器中的失效特征点。f_manager: feature_manager
                ROS_INFO("Initialization finish!");
                last_R = Rs[WINDOW_SIZE]; // 更新上一帧的旋转和平移。
                last_P = Ps[WINDOW_SIZE];
                last_R0 = Rs[0];
                last_P0 = Ps[0];    
            }
            else
            {
                slideWindow();
            }
        }
        else
            frame_count++;
    }
    else // 执行非初始化的一系列操作
    {
        TicToc t_solve;
        solveOdometry();
        ROS_DEBUG("solver costs: %fms", t_solve.toc());

        if (failureDetection()) // 失效检测
        {
            ROS_WARN("failure detection!");
            failure_occur = 1;
            clearState();
            setParameter();
            ROS_WARN("system reboot!");
            return;
        }

        TicToc t_margin;
        slideWindow();
        f_manager.removeFailures();
        ROS_DEBUG("marginalization costs: %fms", t_margin.toc());
        key_poses.clear();
        for (int i = 0; i <= WINDOW_SIZE; i++)
            key_poses.push_back(Ps[i]);

        last_R = Rs[WINDOW_SIZE];
        last_P = Ps[WINDOW_SIZE];
        last_R0 = Rs[0];
        last_P0 = Ps[0];
    }
}

void Estimator::inputEphem(EphemBasePtr ephem_ptr)
{
    double toe = time2sec(ephem_ptr->toe); // 时间转换为秒
    // if a new ephemeris comes
    if (sat2time_index.count(ephem_ptr->sat) == 0 || sat2time_index.at(ephem_ptr->sat).count(toe) == 0)
    {
        sat2ephem[ephem_ptr->sat].emplace_back(ephem_ptr);
        sat2time_index[ephem_ptr->sat].emplace(toe, sat2ephem.at(ephem_ptr->sat).size()-1);
    }
}

void Estimator::inputIonoParams(double ts, const std::vector<double> &iono_params)
{
    if (iono_params.size() != 8)    return;

    // update ionosphere parameters
    latest_gnss_iono_params.clear();
    std::copy(iono_params.begin(), iono_params.end(), std::back_inserter(latest_gnss_iono_params)); // 不断更新电离层参数
}

void Estimator::inputGNSSTimeDiff(const double t_diff)
{
    diff_t_gnss_local = t_diff;
}

// 【GNSS信息滤波器】
// 根据各种条件，获取可用的GNSS观测值数据
void Estimator::processGNSS(const std::vector<ObsPtr> &gnss_meas)
{
    std::vector<ObsPtr> valid_meas; // 有用的测量数据
    std::vector<EphemBasePtr> valid_ephems; // 有用的星历数据

    for (auto obs : gnss_meas) // 遍历所有的GNSS测量数据，得到有用信息
    {
        // Step 1：根据四个导航系统特征，确定当前导航系统
        // filter according to system
        uint32_t sys = satsys(obs->sat, NULL);
        if (sys != SYS_GPS && sys != SYS_GLO && sys != SYS_GAL && sys != SYS_BDS)
            continue;

        // if not got cooresponding ephemeris yet
        // Step 2：确定当前时刻观测到的卫星数量
        if (sat2ephem.count(obs->sat) == 0)
            continue;
        
        // Step 3：如果观测数据中频率值不为空，选择L1频率，并确定相应的观测值和星历信息
        if (obs->freqs.empty())    continue;       // no valid signal measurement 
        int freq_idx = -1;
        L1_freq(obs, &freq_idx); // 通过 L1_freq 函数获取L1频率的索引
        if (freq_idx < 0)   continue;              // no L1 observation
        
        // Step 4：根据星历的参数时间toe，判断当前星历参数是否有效
        // 一般认为当前的GPS时间在星历参考时间前后的2小时之内（即7200秒之内），这套星历参数被认为是有效的
        // PS：GPS的星历在偶数整点更新（2小时更新一次）
        double obs_time = time2sec(obs->time); // 将观测时间转换为秒
        std::map<double, size_t> time2index = sat2time_index.at(obs->sat); // 获取当前卫星的时间戳和索引的映射。
        double ephem_time = EPH_VALID_SECONDS; // 星历时间
        size_t ephem_index = -1; // 星历索引
        for (auto ti : time2index) // 遍历卫星的时间戳和索引的映射。
        {
            if (std::abs(ti.first - obs_time) < ephem_time) // 如果当前观测时间与星历时间的差值小于 EPH_VALID_SECONDS，更新星历时间和星历索引。
            {
                ephem_time = std::abs(ti.first - obs_time);
                ephem_index = ti.second;
            }
        }
        if (ephem_time >= EPH_VALID_SECONDS) // 超过7200秒，丢弃（星历不可用） 若星历时间超过了 EPH_VALID_SECONDS，则认为星历不再有效，跳过该次循环。
        {
            ROS_WARN("ephemeris not valid anymore");
            // cerr << "ephemeris not valid anymore\n";
            continue;
        }
        const EphemBasePtr &best_ephem = sat2ephem.at(obs->sat).at(ephem_index); // 获取当前卫星对应的最佳星历数据。

        // filter by tracking status
        // Step 5：根据卫星的跟踪状态和伪距、多普勒测量的标准差，判断该测量数据的合法性
        LOG_IF(FATAL, freq_idx < 0) << "No L1 observation found.\n";
        // psr_std：伪距标准差
        // dopp_std：多普勒标准差
        // 如果伪距和多普勒测量中任一数据的标准差超过了规定数值，认为当前测量值不OK，并将当前测量值下的卫星被跟踪（被观测）次数状态归零
        // PSR_STD 其实代表了伪距的不确定度； DOPP_STD 代表了多普勒的不确定度
        if (obs->psr_std[freq_idx]  > GNSS_PSR_STD_THRES ||
            obs->dopp_std[freq_idx] > GNSS_DOPP_STD_THRES) // 判断伪距和多普勒测量的标准差是否超过阈值。如果超过阈值，将该卫星的跟踪状态清零。
        {
            sat_track_status[obs->sat] = 0;
            continue;
        }
        else // 如果伪距和多普勒测量的标准差在阈值内，更新该卫星的跟踪次数。
        {
            if (sat_track_status.count(obs->sat) == 0)
                sat_track_status[obs->sat] = 0;
            ++ sat_track_status[obs->sat]; // 跟踪（观测）次数+1
        }
        // 如果卫星的跟踪次数不足 GNSS_TRACK_NUM_THRES，表示卫星锁定时间不够，跳过该次循环。
        if (sat_track_status[obs->sat] < GNSS_TRACK_NUM_THRES)
            continue;           // not being tracked for enough epochs

        // Step 6：根据观测卫星的仰角
        // filter by elevation angle
        if (gnss_ready) // 如果GNSS已经准备好（当GNSS-VI-Align()对齐之后，gnss_ready为true）
        {
            // 根据卫星的导航系统，获取卫星的ECEF坐标。
            Eigen::Vector3d sat_ecef;
            if (sys == SYS_GLO)
                sat_ecef = geph2pos(obs->time, std::dynamic_pointer_cast<GloEphem>(best_ephem), NULL);
            else
                sat_ecef = eph2pos(obs->time, std::dynamic_pointer_cast<Ephem>(best_ephem), NULL);
            
            // 计算卫星在接收机位置的仰角。这里 ecef_pos 表示接收机的ECEF坐标。
            double azel[2] = {0, M_PI/2.0};
            sat_azel(ecef_pos, sat_ecef, azel); // 通过ENU坐标系下，接收机位置和卫星位置的向量差，可求出仰角
            if (azel[1] < GNSS_ELEVATION_THRES*M_PI/180.0) // GNSS_ELEVATION_THRES=30度，要求仰角不能小于30度
                continue;
        }
        valid_meas.push_back(obs);
        valid_ephems.push_back(best_ephem);
    }
    
    // Step 7：将条件较好的卫星观测数据和星历信息放到全局变量中
    gnss_meas_buf[frame_count] = valid_meas; //原始量测
    gnss_ephem_buf[frame_count] = valid_ephems; //星历
}

bool Estimator::initialStructure()
{
    TicToc t_sfm;
    //check imu observibility
    {
        map<double, ImageFrame>::iterator frame_it;
        Vector3d sum_g;
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            sum_g += tmp_g;
        }
        Vector3d aver_g;
        aver_g = sum_g * 1.0 / ((int)all_image_frame.size() - 1);
        double var = 0;
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            var += (tmp_g - aver_g).transpose() * (tmp_g - aver_g);
            //cout << "frame g " << tmp_g.transpose() << endl;
        }
        var = sqrt(var / ((int)all_image_frame.size() - 1));
        //ROS_WARN("IMU variation %f!", var);
        if(var < 0.25)
        {
            ROS_INFO("IMU excitation not enough!");
            //return false;
        }
    }
    // global sfm
    Quaterniond Q[frame_count + 1];
    Vector3d T[frame_count + 1];
    map<int, Vector3d> sfm_tracked_points;
    vector<SFMFeature> sfm_f;
    for (auto &it_per_id : f_manager.feature)
    {
        int imu_j = it_per_id.start_frame - 1;
        SFMFeature tmp_feature;
        tmp_feature.state = false;
        tmp_feature.id = it_per_id.feature_id;
        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            Vector3d pts_j = it_per_frame.point;
            tmp_feature.observation.push_back(make_pair(imu_j, Eigen::Vector2d{pts_j.x(), pts_j.y()}));
        }
        sfm_f.push_back(tmp_feature);
    } 
    Matrix3d relative_R;
    Vector3d relative_T;
    int l;
    if (!relativePose(relative_R, relative_T, l))
    {
        ROS_INFO("Not enough features or parallax; Move device around");
        return false;
    }
    GlobalSFM sfm;
    if(!sfm.construct(frame_count + 1, Q, T, l,
              relative_R, relative_T,
              sfm_f, sfm_tracked_points))
    {
        ROS_DEBUG("global SFM failed!");
        marginalization_flag = MARGIN_OLD;
        return false;
    }

    //solve pnp for all frame
    map<double, ImageFrame>::iterator frame_it;
    map<int, Vector3d>::iterator it;
    frame_it = all_image_frame.begin( );
    for (int i = 0; frame_it != all_image_frame.end( ); frame_it++)
    {
        // provide initial guess
        cv::Mat r, rvec, t, D, tmp_r;
        if((frame_it->first) == Headers[i].stamp.toSec())
        {
            frame_it->second.is_key_frame = true;
            frame_it->second.R = Q[i].toRotationMatrix() * RIC[0].transpose();
            frame_it->second.T = T[i];
            i++;
            continue;
        }
        if((frame_it->first) > Headers[i].stamp.toSec())
        {
            i++;
        }
        Matrix3d R_inital = (Q[i].inverse()).toRotationMatrix();
        Vector3d P_inital = - R_inital * T[i];
        cv::eigen2cv(R_inital, tmp_r);
        cv::Rodrigues(tmp_r, rvec);
        cv::eigen2cv(P_inital, t);

        frame_it->second.is_key_frame = false;
        vector<cv::Point3f> pts_3_vector;
        vector<cv::Point2f> pts_2_vector;
        for (auto &id_pts : frame_it->second.points)
        {
            int feature_id = id_pts.first;
            for (auto &i_p : id_pts.second)
            {
                it = sfm_tracked_points.find(feature_id);
                if(it != sfm_tracked_points.end())
                {
                    Vector3d world_pts = it->second;
                    cv::Point3f pts_3(world_pts(0), world_pts(1), world_pts(2));
                    pts_3_vector.push_back(pts_3);
                    Vector2d img_pts = i_p.second.head<2>();
                    cv::Point2f pts_2(img_pts(0), img_pts(1));
                    pts_2_vector.push_back(pts_2);
                }
            }
        }
        cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);     
        if(pts_3_vector.size() < 6)
        {
            cout << "pts_3_vector size " << pts_3_vector.size() << endl;
            ROS_DEBUG("Not enough points for solve pnp !");
            return false;
        }
        if (! cv::solvePnP(pts_3_vector, pts_2_vector, K, D, rvec, t, 1))
        {
            ROS_DEBUG("solve pnp fail!");
            return false;
        }
        cv::Rodrigues(rvec, r);
        MatrixXd R_pnp,tmp_R_pnp;
        cv::cv2eigen(r, tmp_R_pnp);
        R_pnp = tmp_R_pnp.transpose();
        MatrixXd T_pnp;
        cv::cv2eigen(t, T_pnp);
        T_pnp = R_pnp * (-T_pnp);
        frame_it->second.R = R_pnp * RIC[0].transpose();
        frame_it->second.T = T_pnp;
    }

    if (!visualInitialAlign())
    {
        ROS_WARN("misalign visual structure with IMU");
        return false;
    }
    return true;
}

bool Estimator::visualInitialAlign()
{
    TicToc t_g;
    VectorXd x;
    //solve scale
    bool result = VisualIMUAlignment(all_image_frame, Bgs, g, x);
    if(!result)
    {
        ROS_DEBUG("solve g failed!");
        return false;
    }

    // change state
    for (int i = 0; i <= frame_count; i++)
    {
        Matrix3d Ri = all_image_frame[Headers[i].stamp.toSec()].R;
        Vector3d Pi = all_image_frame[Headers[i].stamp.toSec()].T;
        Ps[i] = Pi;
        Rs[i] = Ri;
        all_image_frame[Headers[i].stamp.toSec()].is_key_frame = true;
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < dep.size(); i++)
        dep[i] = -1;
    f_manager.clearDepth(dep);

    //triangulate on cam pose , no tic
    Vector3d TIC_TMP[NUM_OF_CAM];
    for(int i = 0; i < NUM_OF_CAM; i++)
        TIC_TMP[i].setZero();
    ric[0] = RIC[0];
    f_manager.setRic(ric);
    f_manager.triangulate(Ps, &(TIC_TMP[0]), &(RIC[0]));

    double s = (x.tail<1>())(0);
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        pre_integrations[i]->repropagate(Vector3d::Zero(), Bgs[i]);
    }
    for (int i = frame_count; i >= 0; i--)
        Ps[i] = s * Ps[i] - Rs[i] * TIC[0] - (s * Ps[0] - Rs[0] * TIC[0]);
    int kv = -1;
    map<double, ImageFrame>::iterator frame_i;
    for (frame_i = all_image_frame.begin(); frame_i != all_image_frame.end(); frame_i++)
    {
        if(frame_i->second.is_key_frame)
        {
            kv++;
            Vs[kv] = frame_i->second.R * x.segment<3>(kv * 3);
        }
    }
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        it_per_id.estimated_depth *= s;
    }

    Matrix3d R0 = Utility::g2R(g);
    double yaw = Utility::R2ypr(R0 * Rs[0]).x();
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
    g = R0 * g;
    //Matrix3d rot_diff = R0 * Rs[0].transpose();
    Matrix3d rot_diff = R0;
    for (int i = 0; i <= frame_count; i++)
    {
        Ps[i] = rot_diff * Ps[i];
        Rs[i] = rot_diff * Rs[i];
        Vs[i] = rot_diff * Vs[i];
    }

    ROS_DEBUG_STREAM("g0     " << g.transpose());
    ROS_DEBUG_STREAM("my R0  " << Utility::R2ypr(Rs[0]).transpose());

    return true;
}

/**
 * @brief GNSS和VIO进行联合优化，在VIO初始化完成之后
 * **/
bool Estimator::GNSSVIAlign()
{
    static int locked_sat_num = 0;
    
    // Step 1：VIO初始化未完成，不进行GNSS和VIO的对齐
    if (solver_flag == INITIAL)     // visual-inertial not initialized
        return false;
    
    // 如果已经完成GNSS和VIO的对齐，直接退出返回true
    if (gnss_ready)                 // GNSS-VI already initialized
        return true;
    
     // 如果GNSS的观测数据为空或者每帧对应下的观测数据数量太少（小于10），直接退出返回false，不继续进行对齐
    // std::cout << "----------------------------" << std::endl;
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
    {
        if (gnss_meas_buf[i].empty() || gnss_meas_buf[i].size() < 10){ // 这个条件有些苛刻，可以考虑将10改小（eg. 4）,或者直接去掉 gnss_meas_buf[i].size() < 10 这个条件
            if(static_cast<int>(gnss_meas_buf[i].size()) != locked_sat_num) {
                locked_sat_num = gnss_meas_buf[i].size();
            }
            return false;
        } else {
            if(static_cast<int>(gnss_meas_buf[i].size()) != locked_sat_num) {
                locked_sat_num = gnss_meas_buf[i].size();
            }
        }
    }

    // Step 2：检测水平方向的速度激励,先计算一个滑窗内水平速度的平均值，并通过该平均值的模（magnitude）是否小于0.3来判断速度激励是否足够
    // check horizontal velocity excitation
    Eigen::Vector2d avg_hor_vel(0.0, 0.0); // 用于存储水平方向上的速度平均值，初始值设为 (0.0, 0.0)。
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
        avg_hor_vel += Vs[i].head<2>().cwiseAbs(); // 将当前帧的水平速度的绝对值（只考虑方向）累加到 avg_hor_vel 上。 Vs[i]是第i帧的速度向量；head<2>() 提取速度向量的前两个元素，即水平方向上的速度。 cwiseAbs() 对速度向量的每个元素取绝对值。
    avg_hor_vel /= (WINDOW_SIZE+1); //当前滑窗中， 水平方向上的速度平均值
    if (avg_hor_vel.norm() < 0.3)  // 如果平均速度的模小于0.3，认为速度激励不够
    {
        // 在GNSS论文中有提到，将行人正常移动的平均速度认为是0.3m/s
        // 应该是从这个数值入手，判断当前手持设备是否正常移动
        // TODO：速度激励不够，即移动过慢，会造成ECEF系到ENU系变换没有意义？？？
        std::cerr << "velocity excitation not enough for GNSS-VI alignment.\n";
        return false;
    }

    // 获取当前窗口内的GNSS观测数据和星历数据
    std::vector<std::vector<ObsPtr>> curr_gnss_meas_buf; // 当前GNSS观测值
    std::vector<std::vector<EphemBasePtr>> curr_gnss_ephem_buf; // 当前星历信息
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i) 
    {
        curr_gnss_meas_buf.push_back(gnss_meas_buf[i]);
        curr_gnss_ephem_buf.push_back(gnss_ephem_buf[i]);
    }

    // Step 3：进行GNSS和VIO的对齐
    GNSSVIInitializer gnss_vi_initializer(curr_gnss_meas_buf, curr_gnss_ephem_buf, latest_gnss_iono_params); // latest_gnss_iono_params是电离层参数

    // 1. get a rough global location
    // Step 3.1：基于SPP算法，得到接收器的一个在ECEF系下的粗糙伪距位置和四个导航系统中卫星的时钟钟差
    Eigen::Matrix<double, 7, 1> rough_xyzt;
    rough_xyzt.setZero();
    if (!gnss_vi_initializer.coarse_localization(rough_xyzt))
    {
        std::cerr << "Fail to obtain a coarse location.\n";
        return false;
    }

    // 2. perform yaw alignment
    // Step 3.2：local world frame和ENU系的对齐，即yaw的求解，毕竟这两个坐标系的Z轴是完完全全重合的
    // 使用多普勒测量，比伪距测量精度高一个数量级
    std::vector<Eigen::Vector3d> local_vs; // local world frame中的速度。计算得是从body 相对于 local world 的速度,其实就是VIO的速度
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
        local_vs.push_back(Vs[i]);
    Eigen::Vector3d rough_anchor_ecef = rough_xyzt.head<3>(); // 接收机在ECEF下的粗糙位置
    double aligned_yaw = 0; // 待校正的yaw角
    double aligned_rcv_ddt = 0; // 接收机时钟钟差变化率
    if (!gnss_vi_initializer.yaw_alignment(local_vs, rough_anchor_ecef, aligned_yaw, aligned_rcv_ddt)) // rough_anchor_ecef ： 接收机的粗位置； aligned_yaw：偏航角
    {
        std::cerr << "Fail to align ENU and local frames.\n";
        return false;
    }
    // std::cout << "aligned_yaw is " << aligned_yaw*180.0/M_PI << '\n';

    // 3. perform anchor refinement
    // Step 3.2：瞄点位置重优化
    std::vector<Eigen::Vector3d> local_ps; // VIO在local world frame的位置
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
        local_ps.push_back(Ps[i]);
    Eigen::Matrix<double, 7, 1> refined_xyzt; // anchor point位置和四个导航的卫星钟差
    refined_xyzt.setZero();
    if (!gnss_vi_initializer.anchor_refinement(local_ps, aligned_yaw, 
        aligned_rcv_ddt, rough_xyzt, refined_xyzt))
    {
        std::cerr << "Fail to refine anchor point.\n";
        return false;
    }
    // std::cout << "refined anchor point is " << std::setprecision(20) 
    //           << refined_xyzt.head<3>().transpose() << '\n';

    // restore GNSS states
    uint32_t one_observed_sys = static_cast<uint32_t>(-1); // 被观测到的导航系统
    for (uint32_t k = 0; k < 4; ++k)
    {
        if (rough_xyzt(k+3) != 0)
        {
            one_observed_sys = k;
            break;
        }
    }
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
    {
        para_rcv_ddt[i] = aligned_rcv_ddt;
        for (uint32_t k = 0; k < 4; ++k)
        {
            if (rough_xyzt(k+3) == 0)
                para_rcv_dt[i*4+k] = refined_xyzt(3+one_observed_sys) + aligned_rcv_ddt * i;
            else
                para_rcv_dt[i*4+k] = refined_xyzt(3+k) + aligned_rcv_ddt * i;
        }
    }
    anc_ecef = refined_xyzt.head<3>();
    R_ecef_enu = ecef2rotation(anc_ecef);

    yaw_enu_local = aligned_yaw;

    return true;
}

void Estimator::updateGNSSStatistics()
{
    R_enu_local = Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ());
    enu_pos = R_enu_local * Ps[WINDOW_SIZE];
    enu_vel = R_enu_local * Vs[WINDOW_SIZE];
    enu_ypr = Utility::R2ypr(R_enu_local*Rs[WINDOW_SIZE]);
    ecef_pos = anc_ecef + R_ecef_enu * enu_pos;
}


bool Estimator::relativePose(Matrix3d &relative_R, Vector3d &relative_T, int &l)
{
    // find previous frame which contians enough correspondance and parallex with newest frame
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        vector<pair<Vector3d, Vector3d>> corres;
        corres = f_manager.getCorresponding(i, WINDOW_SIZE);
        if (corres.size() > 20)
        {
            double sum_parallax = 0;
            double average_parallax;
            for (int j = 0; j < int(corres.size()); j++)
            {
                Vector2d pts_0(corres[j].first(0), corres[j].first(1));
                Vector2d pts_1(corres[j].second(0), corres[j].second(1));
                double parallax = (pts_0 - pts_1).norm();
                sum_parallax = sum_parallax + parallax;

            }
            average_parallax = 1.0 * sum_parallax / int(corres.size());
            if(average_parallax * 460 > 30 && m_estimator.solveRelativeRT(corres, relative_R, relative_T))
            {
                l = i;
                ROS_DEBUG("average_parallax %f choose l %d and newest frame to triangulate the whole structure", average_parallax * 460, l);
                return true;
            }
        }
    }
    return false;
}

void Estimator::solveOdometry()
{
    if (frame_count < WINDOW_SIZE)
        return;
    if (solver_flag == NON_LINEAR)
    {
        TicToc t_tri;
        f_manager.triangulate(Ps, tic, ric);
        ROS_DEBUG("triangulation costs %f", t_tri.toc());
        optimization();
        if (GNSS_ENABLE.load())
        {
            if (!gnss_ready)
            {
                gnss_ready = GNSSVIAlign();
                
            }
            if (gnss_ready)
            {
                updateGNSSStatistics();
            }
        }
    }
}

void Estimator::vector2double()
{
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        para_Pose[i][0] = Ps[i].x();
        para_Pose[i][1] = Ps[i].y();
        para_Pose[i][2] = Ps[i].z();
        Quaterniond q{Rs[i]};
        para_Pose[i][3] = q.x();
        para_Pose[i][4] = q.y();
        para_Pose[i][5] = q.z();
        para_Pose[i][6] = q.w();

        para_SpeedBias[i][0] = Vs[i].x();
        para_SpeedBias[i][1] = Vs[i].y();
        para_SpeedBias[i][2] = Vs[i].z();

        para_SpeedBias[i][3] = Bas[i].x();
        para_SpeedBias[i][4] = Bas[i].y();
        para_SpeedBias[i][5] = Bas[i].z();

        para_SpeedBias[i][6] = Bgs[i].x();
        para_SpeedBias[i][7] = Bgs[i].y();
        para_SpeedBias[i][8] = Bgs[i].z();
    }
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        para_Ex_Pose[i][0] = tic[i].x();
        para_Ex_Pose[i][1] = tic[i].y();
        para_Ex_Pose[i][2] = tic[i].z();
        Quaterniond q{ric[i]};
        para_Ex_Pose[i][3] = q.x();
        para_Ex_Pose[i][4] = q.y();
        para_Ex_Pose[i][5] = q.z();
        para_Ex_Pose[i][6] = q.w();
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        para_Feature[i][0] = dep(i);

    if (ESTIMATE_TD)
        para_Td[0][0] = td;
    
    para_yaw_enu_local[0] = yaw_enu_local;
    for (uint32_t k = 0; k < 3; ++k)
        para_anc_ecef[k] = anc_ecef(k);
}

void Estimator::double2vector()
{
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {

        Rs[i] = Quaterniond(para_Pose[i][6], para_Pose[i][3], 
                            para_Pose[i][4], para_Pose[i][5]).normalized().toRotationMatrix();
        
        Ps[i] = Vector3d(para_Pose[i][0], para_Pose[i][1], para_Pose[i][2]);

        Vs[i] = Vector3d(para_SpeedBias[i][0], para_SpeedBias[i][1], para_SpeedBias[i][2]);

        Bas[i] = Vector3d(para_SpeedBias[i][3], para_SpeedBias[i][4], para_SpeedBias[i][5]);

        Bgs[i] = Vector3d(para_SpeedBias[i][6], para_SpeedBias[i][7], para_SpeedBias[i][8]);
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d(para_Ex_Pose[i][0], para_Ex_Pose[i][1], para_Ex_Pose[i][2]);
        ric[i] = Quaterniond(para_Ex_Pose[i][6], para_Ex_Pose[i][3],
                             para_Ex_Pose[i][4], para_Ex_Pose[i][5]).normalized().toRotationMatrix();
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        dep(i) = para_Feature[i][0];
    f_manager.setDepth(dep);
    if (ESTIMATE_TD)
        td = para_Td[0][0];
    
    if (gnss_ready)
    {
        yaw_enu_local = para_yaw_enu_local[0];
        for (uint32_t k = 0; k < 3; ++k)
            anc_ecef(k) = para_anc_ecef[k];
        R_ecef_enu = ecef2rotation(anc_ecef);
    }
}

bool Estimator::failureDetection()
{
    if (f_manager.last_track_num < 2)
    {
        ROS_INFO(" little feature %d", f_manager.last_track_num);
        //return true;
    }
    if (Bas[WINDOW_SIZE].norm() > 5.0) // [Damon correct 5.0] 原来是2.5
    {
        // ROS_INFO(" big IMU acc bias estimation %f", Bas[WINDOW_SIZE].norm());
        // return true;
    }
    if (Bgs[WINDOW_SIZE].norm() > 1.0)
    {
        ROS_INFO(" big IMU gyr bias estimation %f", Bgs[WINDOW_SIZE].norm());
        // return true;
    }
    /*
    if (tic(0) > 1)
    {
        ROS_INFO(" big extri param estimation %d", tic(0) > 1);
        return true;
    }
    */
    // Vector3d tmp_P = Ps[WINDOW_SIZE];
    // if ((tmp_P - last_P).norm() > 5) 
    // {
    //     ROS_INFO(" big translation");
    //     return true;
    // }
    // if (abs(tmp_P.z() - last_P.z()) > 1)
    // {
    //     ROS_INFO(" big z translation");
    //     return true; 
    // }
    // Matrix3d tmp_R = Rs[WINDOW_SIZE];
    // Matrix3d delta_R = tmp_R.transpose() * last_R;
    // Quaterniond delta_Q(delta_R);
    // double delta_angle;
    // delta_angle = acos(delta_Q.w()) * 2.0 / 3.14 * 180.0;
    // if (delta_angle > 50)
    // {
    //     ROS_INFO(" big delta_angle ");
    //     //return true;
    // }
    return false;
}

void Estimator::optimization()
{
    ceres::Problem problem; // 声明了一个 ceres::Problem 类型的对象 problem，用于表示优化问题
    ceres::LossFunction *loss_function; // 声明了一个指向 ceres::LossFunction 类型的指针 loss_function
    //loss_function = new ceres::HuberLoss(1.0);
    loss_function = new ceres::CauchyLoss(1.0);
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization(); // 优化中的位姿（Pose）参数块
        problem.AddParameterBlock(para_Pose[i], SIZE_POSE, local_parameterization); // 位置3维 四元数4维
        problem.AddParameterBlock(para_SpeedBias[i], SIZE_SPEEDBIAS); // 速度3维 和 IMU零偏 ba 3维  bg3维，共9维
    }
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Ex_Pose[i], SIZE_POSE, local_parameterization);
        if (!ESTIMATE_EXTRINSIC)
        {
            ROS_DEBUG("fix extinsic param");
            problem.SetParameterBlockConstant(para_Ex_Pose[i]); // 将优化问题中的某个参数块设置为常数，表示在优化过程中该参数块的值不会被更新
        }
        else
            ROS_DEBUG("estimate extinsic param");
    }
    if (ESTIMATE_TD)
    {
        problem.AddParameterBlock(para_Td[0], 1);
    }

    if (gnss_ready)
    {
        problem.AddParameterBlock(para_yaw_enu_local, 1);
        Eigen::Vector2d avg_hor_vel(0.0, 0.0); // 水平速度的平均值
        for (uint32_t i = 0; i <= WINDOW_SIZE; ++i)
            avg_hor_vel += Vs[i].head<2>().cwiseAbs(); // head<2>取Vs[i]向量的前两个分量，cwiseAbs()为取绝对值
        avg_hor_vel /= (WINDOW_SIZE+1); //窗口内所有的速度水平分量的平均值
        // cerr << "avg_hor_vel is " << avg_vel << endl;
        if (avg_hor_vel.norm() < 0.3)
        {
            // std::cerr << "velocity excitation not enough, fix yaw angle.\n";
            // 如果水平速度小于 0.3，意味着当前的速度变化不够强烈，可能没有足够的运动数据来有效估计yaw_enu_local，在这种情况下，不再对其进行优化。
            problem.SetParameterBlockConstant(para_yaw_enu_local); // 
        }

        for (uint32_t i = 0; i <= WINDOW_SIZE; ++i)
        {
            if (gnss_meas_buf[i].size() < 10) // 该时刻的 GNSS 测量数据数量小于 10， 不再优化 para_yaw_enu_local。
                problem.SetParameterBlockConstant(para_yaw_enu_local);
        }
        
        problem.AddParameterBlock(para_anc_ecef, 3); //anchor的ECEF坐标
        // problem.SetParameterBlockConstant(para_anc_ecef);

        for (uint32_t i = 0; i <= WINDOW_SIZE; ++i) //接收机钟差
        {
            for (uint32_t k = 0; k < 4; ++k)
                problem.AddParameterBlock(para_rcv_dt+i*4+k, 1);
            problem.AddParameterBlock(para_rcv_ddt+i, 1);
        }
    }

    TicToc t_whole, t_prepare;
    vector2double(); //将estimator的值传入估计器临时变量

    if (first_optimization) //首次优化时传入 pose_anchor 因子
    {
        std::vector<double> anchor_value;
        for (uint32_t k = 0; k < 7; ++k)
            anchor_value.push_back(para_Pose[0][k]);
        PoseAnchorFactor *pose_anchor_factor = new PoseAnchorFactor(anchor_value);
        problem.AddResidualBlock(pose_anchor_factor, NULL, para_Pose[0]);
        first_optimization = false;
    }

    if (last_marginalization_info) // 边缘化残差块 last_marginalization_info是个指针，判空指针
    {
        // construct new marginlization_factor
        MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
        problem.AddResidualBlock(marginalization_factor, NULL,
                                 last_marginalization_parameter_blocks);
    }

    for (int i = 0; i < WINDOW_SIZE; i++) // 预积分残差块
    {
        int j = i + 1;
        if (pre_integrations[j]->sum_dt > 10.0)
            continue;
        IMUFactor* imu_factor = new IMUFactor(pre_integrations[j]);
        problem.AddResidualBlock(imu_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j]);
    }

    if (gnss_ready)
    {
        for(int i = 0; i <= WINDOW_SIZE; ++i)
        {
            // cerr << "size of gnss_meas_buf[" << i << "] is " << gnss_meas_buf[i].size() << endl;
            const std::vector<ObsPtr> &curr_obs = gnss_meas_buf[i];
            const std::vector<EphemBasePtr> &curr_ephem = gnss_ephem_buf[i];

            for (uint32_t j = 0; j < curr_obs.size(); ++j) // 伪距 多普勒残差 
            {
                const uint32_t sys = satsys(curr_obs[j]->sat, NULL);
                const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);

                int lower_idx = -1;
                const double obs_local_ts = time2sec(curr_obs[j]->time) - diff_t_gnss_local;
                if (Headers[i].stamp.toSec() > obs_local_ts)
                    lower_idx = (i==0? 0 : i-1);
                else
                    lower_idx = (i==WINDOW_SIZE? WINDOW_SIZE-1 : i);
                const double lower_ts = Headers[lower_idx].stamp.toSec();
                const double upper_ts = Headers[lower_idx+1].stamp.toSec();

                const double ts_ratio = (upper_ts-obs_local_ts) / (upper_ts-lower_ts);
                GnssPsrDoppFactor *gnss_factor = new GnssPsrDoppFactor(curr_obs[j], 
                    curr_ephem[j], latest_gnss_iono_params, ts_ratio);
                problem.AddResidualBlock(gnss_factor, NULL, para_Pose[lower_idx], 
                    para_SpeedBias[lower_idx], para_Pose[lower_idx+1], para_SpeedBias[lower_idx+1],
                    para_rcv_dt+i*4+sys_idx, para_rcv_ddt+i, para_yaw_enu_local, para_anc_ecef);
            }
        }

        // build relationship between rcv_dt and rcv_ddt
        for (size_t k = 0; k < 4; ++k)
        {
            for (uint32_t i = 0; i < WINDOW_SIZE; ++i)
            {
                const double gnss_dt = Headers[i+1].stamp.toSec() - Headers[i].stamp.toSec();
                DtDdtFactor *dt_ddt_factor = new DtDdtFactor(gnss_dt);
                problem.AddResidualBlock(dt_ddt_factor, NULL, para_rcv_dt+i*4+k, 
                    para_rcv_dt+(i+1)*4+k, para_rcv_ddt+i, para_rcv_ddt+i+1);
            }
        }

        // add rcv_ddt smooth factor
        for (int i = 0; i < WINDOW_SIZE; ++i)
        {
            DdtSmoothFactor *ddt_smooth_factor = new DdtSmoothFactor(GNSS_DDT_WEIGHT);
            problem.AddResidualBlock(ddt_smooth_factor, NULL, para_rcv_ddt+i, para_rcv_ddt+i+1);
        }
    }

    int f_m_cnt = 0;
    int feature_index = -1;
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
 
        ++feature_index;

        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
        
        Vector3d pts_i = it_per_id.feature_per_frame[0].point;

        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            if (imu_i == imu_j)
            {
                continue;
            }
            Vector3d pts_j = it_per_frame.point;
            if (ESTIMATE_TD)
            {
                    ProjectionTdFactor *f_td = new ProjectionTdFactor(pts_i, pts_j, 
                        it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                        it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                    problem.AddResidualBlock(f_td, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]);
            }
            else
            {
                ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                problem.AddResidualBlock(f, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index]);
            }
            f_m_cnt++;
        }
    }

    ROS_DEBUG("visual measurement count: %d", f_m_cnt);
    ROS_DEBUG("prepare for ceres: %f", t_prepare.toc());

    ceres::Solver::Options options; // 配置类，用于设置求解器的各种选项
    options.linear_solver_type = ceres::DENSE_SCHUR; // 设置线性求解器类型为DENSE_SCHUR，它是基于Schur补的密集矩阵线性求解方法，适用于相对较小的优化问题。通过Schur分解求解线性系统
    //options.num_threads = 2;
    options.trust_region_strategy_type = ceres::DOGLEG; //  设置信赖域策略为 DOGLEG，这是一种常见的信赖域方法，通过在信赖域内选择步长来保证优化的稳定性。
    options.max_num_iterations = NUM_ITERATIONS; // 设置最大迭代次数，NUM_ITERATIONS 应该是事先定义的常量，表示优化的最大迭代次数。default=8
    //options.use_explicit_schur_complement = true;
    // options.minimizer_progress_to_stdout = true;
    options.use_nonmonotonic_steps = true; // 允许使用非单调步长，即优化过程中的步长可以不一定是单调减小的。这样可以增加优化的灵活性，避免陷入局部最小值。
    if (marginalization_flag == MARGIN_OLD) // 根据 marginalization_flag 的值，设置最大求解时间。
        options.max_solver_time_in_seconds = SOLVER_TIME * 4.0 / 5.0; // 如果 marginalization_flag 等于 MARGIN_OLD，则设置求解器时间为 SOLVER_TIME 的 80%
    else
        options.max_solver_time_in_seconds = SOLVER_TIME; //否则，使用完整的 SOLVER_TIME 时间。
    
    TicToc t_solver; //记录优化过程的时间
    ceres::Solver::Summary summary; // 用于保存优化过程中的总结信息
    ceres::Solve(options, &problem, &summary); //Ceres求解器核心函数，调用这个函数会执行优化并返回结果。
    // cout << summary.BriefReport() << endl;
    // cout << summary.FullReport() << endl;
    ROS_DEBUG("Iterations : %d", static_cast<int>(summary.iterations.size()));
    ROS_DEBUG("solver costs: %f", t_solver.toc());

    while(para_yaw_enu_local[0] > M_PI)   para_yaw_enu_local[0] -= 2.0*M_PI;
    while(para_yaw_enu_local[0] < -M_PI)  para_yaw_enu_local[0] += 2.0*M_PI;
    // std::cout << "yaw is " << para_yaw_enu_local[0]*180/M_PI << std::endl;

    double2vector(); // 将优化器求解结果传入estimator

    TicToc t_whole_marginalization;
    if (marginalization_flag == MARGIN_OLD)
    {
        MarginalizationInfo *marginalization_info = new MarginalizationInfo();
        vector2double();

        if (last_marginalization_info)
        {
            vector<int> drop_set;
            //遍历 last_marginalization_parameter_blocks，确定哪些参数需要被边缘化（通常是历史帧的位姿 para_Pose[0] 和速度偏差 para_SpeedBias[0]）
            for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
            {
                if (last_marginalization_parameter_blocks[i] == para_Pose[0] ||
                    last_marginalization_parameter_blocks[i] == para_SpeedBias[0])
                    drop_set.push_back(i);
            }
            // construct new marginlization_factor
            MarginalizationFactor *marginalization_factor = new MarginalizationFactor(
                last_marginalization_info);
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(
                marginalization_factor, NULL, last_marginalization_parameter_blocks, drop_set);
            marginalization_info->addResidualBlockInfo(residual_block_info);
        }
        else //如果没有找到边缘化信息（即 last_marginalization_info 为 nullptr），则创建一个位置锚定因子 PoseAnchorFactor，并将其添加到 marginalization_info 中。这个因子用来将当前帧的位置作为锚点，保证优化过程中的稳定性。
        {
            std::vector<double> anchor_value;
            for (uint32_t k = 0; k < 7; ++k)
                anchor_value.push_back(para_Pose[0][k]);
            PoseAnchorFactor *pose_anchor_factor = new PoseAnchorFactor(anchor_value);
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(pose_anchor_factor, 
                NULL, vector<double *>{para_Pose[0]}, vector<int>{0});
            marginalization_info->addResidualBlockInfo(residual_block_info);
        }

        {
            if (pre_integrations[1]->sum_dt < 10.0) // 如果 IMU 数据的时间间隔 pre_integrations[1]->sum_dt 小于 10 秒，则添加一个 IMUFactor 来表示 IMU 数据的约束。
            {
                IMUFactor* imu_factor = new IMUFactor(pre_integrations[1]);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(imu_factor, NULL,
                                                                           vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], para_SpeedBias[1]},
                                                                           vector<int>{0, 1});
                marginalization_info->addResidualBlockInfo(residual_block_info);
            }
        }

        if (gnss_ready)
        {
            for (uint32_t j = 0; j < gnss_meas_buf[0].size(); ++j)
            {
                const uint32_t sys = satsys(gnss_meas_buf[0][j]->sat, NULL);
                const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);

                const double obs_local_ts = time2sec(gnss_meas_buf[0][j]->time) - diff_t_gnss_local;
                const double lower_ts = Headers[0].stamp.toSec();
                const double upper_ts = Headers[1].stamp.toSec();
                const double ts_ratio = (upper_ts-obs_local_ts) / (upper_ts-lower_ts);

                GnssPsrDoppFactor *gnss_factor = new GnssPsrDoppFactor(gnss_meas_buf[0][j], 
                    gnss_ephem_buf[0][j], latest_gnss_iono_params, ts_ratio);
                ResidualBlockInfo *psr_dopp_residual_block_info = new ResidualBlockInfo(gnss_factor, NULL,
                    vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], 
                        para_SpeedBias[1],para_rcv_dt+sys_idx, para_rcv_ddt, 
                        para_yaw_enu_local, para_anc_ecef},
                    vector<int>{0, 1, 4, 5});
                marginalization_info->addResidualBlockInfo(psr_dopp_residual_block_info);
            }

            const double gnss_dt = Headers[1].stamp.toSec() - Headers[0].stamp.toSec();
            for (size_t k = 0; k < 4; ++k)
            {
                DtDdtFactor *dt_ddt_factor = new DtDdtFactor(gnss_dt);
                ResidualBlockInfo *dt_ddt_residual_block_info = new ResidualBlockInfo(dt_ddt_factor, NULL,
                    vector<double *>{para_rcv_dt+k, para_rcv_dt+4+k, para_rcv_ddt, para_rcv_ddt+1}, 
                    vector<int>{0, 2});
                marginalization_info->addResidualBlockInfo(dt_ddt_residual_block_info);
            }

            // margin rcv_ddt smooth factor
            DdtSmoothFactor *ddt_smooth_factor = new DdtSmoothFactor(GNSS_DDT_WEIGHT);
            ResidualBlockInfo *ddt_smooth_residual_block_info = new ResidualBlockInfo(ddt_smooth_factor, NULL,
                    vector<double *>{para_rcv_ddt, para_rcv_ddt+1}, vector<int>{0});
            marginalization_info->addResidualBlockInfo(ddt_smooth_residual_block_info);
        }

        { //对于每个特征点，判断它是否在有效的帧范围内（至少有两个观测，并且起始帧小于窗口大小减去2）。如果满足条件，则为每个特征点添加一个 ProjectionFactor 或 ProjectionTdFactor，它们表示特征点在不同帧中的投影关系。
            int feature_index = -1;
            for (auto &it_per_id : f_manager.feature)
            {
                it_per_id.used_num = it_per_id.feature_per_frame.size();
                if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
                    continue;

                ++feature_index;

                int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
                if (imu_i != 0)
                    continue;

                Vector3d pts_i = it_per_id.feature_per_frame[0].point;

                for (auto &it_per_frame : it_per_id.feature_per_frame)
                {
                    imu_j++;
                    if (imu_i == imu_j)
                        continue;

                    Vector3d pts_j = it_per_frame.point;
                    if (ESTIMATE_TD)
                    {
                        ProjectionTdFactor *f_td = new ProjectionTdFactor(pts_i, pts_j, 
                            it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                            it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                        ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f_td, 
                            loss_function, vector<double *>{para_Pose[imu_i], para_Pose[imu_j], 
                                para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]},
                            vector<int>{0, 3});
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                    else
                    {
                        ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                        ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, 
                            loss_function, vector<double *>{para_Pose[imu_i], para_Pose[imu_j], 
                                para_Ex_Pose[0], para_Feature[feature_index]},
                            vector<int>{0, 3});
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                }
            }
        }

        TicToc t_pre_margin;
        marginalization_info->preMarginalize();
        ROS_DEBUG("pre marginalization %f ms", t_pre_margin.toc());
        
        TicToc t_margin;
        marginalization_info->marginalize();
        ROS_DEBUG("marginalization %f ms", t_margin.toc());

        std::unordered_map<long, double *> addr_shift;
        for (int i = 1; i <= WINDOW_SIZE; i++)
        {
            addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
            addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
            for (uint32_t k = 0; k < 4; ++k)
                addr_shift[reinterpret_cast<long>(para_rcv_dt+i*4+k)] = para_rcv_dt+(i-1)*4+k;
            addr_shift[reinterpret_cast<long>(para_rcv_ddt+i)] = para_rcv_ddt+i-1;
        }
        for (int i = 0; i < NUM_OF_CAM; i++)
            addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];
        if (ESTIMATE_TD)
        {
            addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
        }
        addr_shift[reinterpret_cast<long>(para_yaw_enu_local)] = para_yaw_enu_local;
        addr_shift[reinterpret_cast<long>(para_anc_ecef)] = para_anc_ecef;
        vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);

        if (last_marginalization_info)
            delete last_marginalization_info;
        last_marginalization_info = marginalization_info;
        last_marginalization_parameter_blocks = parameter_blocks;
        
    }
    else
    {
        if (last_marginalization_info &&
            std::count(std::begin(last_marginalization_parameter_blocks), std::end(last_marginalization_parameter_blocks), para_Pose[WINDOW_SIZE - 1]))
        {

            MarginalizationInfo *marginalization_info = new MarginalizationInfo();
            vector2double();
            if (last_marginalization_info)
            {
                vector<int> drop_set;
                for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
                {
                    ROS_ASSERT(last_marginalization_parameter_blocks[i] != para_SpeedBias[WINDOW_SIZE - 1]);
                    if (last_marginalization_parameter_blocks[i] == para_Pose[WINDOW_SIZE - 1])
                        drop_set.push_back(i);
                }
                // construct new marginlization_factor
                MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                               last_marginalization_parameter_blocks,
                                                                               drop_set);

                marginalization_info->addResidualBlockInfo(residual_block_info);
            }

            TicToc t_pre_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->preMarginalize();
            ROS_DEBUG("end pre marginalization, %f ms", t_pre_margin.toc());

            TicToc t_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->marginalize();
            ROS_DEBUG("end marginalization, %f ms", t_margin.toc());
            
            std::unordered_map<long, double *> addr_shift;
            for (int i = 0; i <= WINDOW_SIZE; i++)
            {
                if (i == WINDOW_SIZE - 1)
                    continue;
                else if (i == WINDOW_SIZE)
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
                    addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
                    for (uint32_t k = 0; k < 4; ++k)
                        addr_shift[reinterpret_cast<long>(para_rcv_dt+i*4+k)] = para_rcv_dt+(i-1)*4+k;
                    addr_shift[reinterpret_cast<long>(para_rcv_ddt+i)] = para_rcv_ddt+i-1;
                }
                else
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i];
                    addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i];
                    for (uint32_t k = 0; k < 4; ++k)
                        addr_shift[reinterpret_cast<long>(para_rcv_dt+i*4+k)] = para_rcv_dt+i*4+k;
                    addr_shift[reinterpret_cast<long>(para_rcv_ddt+i)] = para_rcv_ddt+i;
                }
            }
            for (int i = 0; i < NUM_OF_CAM; i++)
                addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];
            if (ESTIMATE_TD)
            {
                addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
            }
            addr_shift[reinterpret_cast<long>(para_yaw_enu_local)] = para_yaw_enu_local;
            addr_shift[reinterpret_cast<long>(para_anc_ecef)] = para_anc_ecef;
            vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);
            if (last_marginalization_info)
                delete last_marginalization_info;
            last_marginalization_info = marginalization_info;
            last_marginalization_parameter_blocks = parameter_blocks;
            
        }
    }
    ROS_DEBUG("whole marginalization costs: %f", t_whole_marginalization.toc());
    
    ROS_DEBUG("whole time for ceres: %f", t_whole.toc());
}

void Estimator::slideWindow()
{
    TicToc t_margin;
    if (marginalization_flag == MARGIN_OLD)
    {
        double t_0 = Headers[0].stamp.toSec();
        back_R0 = Rs[0]; //(1) 保存最老帧信息，删除的是滑窗第一帧。
        back_P0 = Ps[0];
        if (frame_count == WINDOW_SIZE)
        {
            //(2) 依次把滑窗内信息前移
            for (int i = 0; i < WINDOW_SIZE; i++)
            {
                Rs[i].swap(Rs[i + 1]);

                std::swap(pre_integrations[i], pre_integrations[i + 1]);

                dt_buf[i].swap(dt_buf[i + 1]);
                linear_acceleration_buf[i].swap(linear_acceleration_buf[i + 1]);
                angular_velocity_buf[i].swap(angular_velocity_buf[i + 1]);

                Headers[i] = Headers[i + 1];
                Ps[i].swap(Ps[i + 1]);
                Vs[i].swap(Vs[i + 1]);
                Bas[i].swap(Bas[i + 1]);
                Bgs[i].swap(Bgs[i + 1]);

                // GNSS related
                gnss_meas_buf[i].swap(gnss_meas_buf[i+1]);
                gnss_ephem_buf[i].swap(gnss_ephem_buf[i+1]);
                for (uint32_t k = 0; k < 4; ++k)
                    para_rcv_dt[i*4+k] = para_rcv_dt[(i+1)*4+k];
                para_rcv_ddt[i] = para_rcv_ddt[i+1];
            }
            //(3) 把滑窗末尾(10帧)信息给最新一帧(11帧),WINDOW_SIZE = 10.
            Headers[WINDOW_SIZE] = Headers[WINDOW_SIZE - 1];
            Ps[WINDOW_SIZE] = Ps[WINDOW_SIZE - 1];
            Vs[WINDOW_SIZE] = Vs[WINDOW_SIZE - 1];
            Rs[WINDOW_SIZE] = Rs[WINDOW_SIZE - 1];
            Bas[WINDOW_SIZE] = Bas[WINDOW_SIZE - 1];
            Bgs[WINDOW_SIZE] = Bgs[WINDOW_SIZE - 1];
            /*注意，在(2)中，已经实现了所有信息的前移，此时，最新一帧已经成为了滑窗中的第10帧，这里只是把原先的最新一帧的信息作为下一次最新一帧的初始值。*/

            // GNSS related
            gnss_meas_buf[WINDOW_SIZE].clear();
            gnss_ephem_buf[WINDOW_SIZE].clear();

            // (4) 新实例化一个IMU预积分对象给下一个最新一帧
            delete pre_integrations[WINDOW_SIZE];
            pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

            //(5) 清空第11帧的buf
            dt_buf[WINDOW_SIZE].clear();
            linear_acceleration_buf[WINDOW_SIZE].clear();
            angular_velocity_buf[WINDOW_SIZE].clear();

            // (6)删除最老帧对应的全部信息.
            if (true || solver_flag == INITIAL)
            {
                map<double, ImageFrame>::iterator it_0;
                it_0 = all_image_frame.find(t_0);
                delete it_0->second.pre_integration;
                it_0->second.pre_integration = nullptr;
 
                for (map<double, ImageFrame>::iterator it = all_image_frame.begin(); it != it_0; ++it)
                {
                    if (it->second.pre_integration)
                        delete it->second.pre_integration;
                    it->second.pre_integration = NULL;
                }

                all_image_frame.erase(all_image_frame.begin(), it_0);
                all_image_frame.erase(t_0);

            }
            slideWindowOld();
        }
    }
    else //(1)取出最新一帧的信息,删除的是滑窗第10帧
    {
        if (frame_count == WINDOW_SIZE)
        {
            for (unsigned int i = 0; i < dt_buf[frame_count].size(); i++)
            {
                double tmp_dt = dt_buf[frame_count][i];
                Vector3d tmp_linear_acceleration = linear_acceleration_buf[frame_count][i];
                Vector3d tmp_angular_velocity = angular_velocity_buf[frame_count][i];

                //(2) 当前帧和前一帧之间的 IMU 预积分转换为当前帧和前二帧之间的 IMU 预积分
                //直接将最后一帧的复制过来不行吗? 注意这里不能复制!!!因为删除了倒数第二帧,所以需要继续预积分到frame_count.
                pre_integrations[frame_count - 1]->push_back(tmp_dt, tmp_linear_acceleration, tmp_angular_velocity);

                dt_buf[frame_count - 1].push_back(tmp_dt);
                linear_acceleration_buf[frame_count - 1].push_back(tmp_linear_acceleration);
                angular_velocity_buf[frame_count - 1].push_back(tmp_angular_velocity);
            }

            //(3) 用最新一帧的信息覆盖上一帧信息
            Headers[frame_count - 1] = Headers[frame_count];
            Ps[frame_count - 1] = Ps[frame_count];
            Vs[frame_count - 1] = Vs[frame_count];
            Rs[frame_count - 1] = Rs[frame_count];
            Bas[frame_count - 1] = Bas[frame_count];
            Bgs[frame_count - 1] = Bgs[frame_count];

            // GNSS related
            gnss_meas_buf[frame_count-1] = gnss_meas_buf[frame_count];
            gnss_ephem_buf[frame_count-1] = gnss_ephem_buf[frame_count];
            for (uint32_t k = 0; k < 4; ++k)
                para_rcv_dt[(frame_count-1)*4+k] = para_rcv_dt[frame_count*4+k];
            para_rcv_ddt[frame_count-1] = para_rcv_ddt[frame_count];
            gnss_meas_buf[frame_count].clear();
            gnss_ephem_buf[frame_count].clear();

            //(4) 因为已经把第11帧的信息覆盖了第10帧，所以现在把第11帧清除
            // 这里delete和new是否会导致计算量大? @todo可以考虑写一个重置的函数
            delete pre_integrations[WINDOW_SIZE];
            pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

            dt_buf[WINDOW_SIZE].clear();
            linear_acceleration_buf[WINDOW_SIZE].clear();
            angular_velocity_buf[WINDOW_SIZE].clear();

            slideWindowNew();
        }
    }
}

// real marginalization is removed in solve_ceres()
void Estimator::slideWindowNew()
{
    sum_of_front++;
    f_manager.removeFront(frame_count);
}
// real marginalization is removed in solve_ceres()
void Estimator::slideWindowOld()
{
    sum_of_back++;

    bool shift_depth = solver_flag == NON_LINEAR ? true : false;
    if (shift_depth)
    {
        Matrix3d R0, R1;
        Vector3d P0, P1;
        R0 = back_R0 * ric[0];
        R1 = Rs[0] * ric[0];
        P0 = back_P0 + back_R0 * tic[0];
        P1 = Ps[0] + Rs[0] * tic[0];
        f_manager.removeBackShiftDepth(R0, P0, R1, P1);
    }
    else
        f_manager.removeBack();
}