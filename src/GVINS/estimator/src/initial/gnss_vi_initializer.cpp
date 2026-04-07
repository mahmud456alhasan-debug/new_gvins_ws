#include "gnss_vi_initializer.h"

GNSSVIInitializer::GNSSVIInitializer(const std::vector<std::vector<ObsPtr>> &gnss_meas_buf_, 
    const std::vector<std::vector<EphemBasePtr>> &gnss_ephem_buf_, const std::vector<double> &iono_params_)
        : gnss_meas_buf(gnss_meas_buf_), gnss_ephem_buf(gnss_ephem_buf_), iono_params(iono_params_)
{
    num_all_meas = 0;
    all_sat_states.clear();
    for (uint32_t i = 0; i < gnss_meas_buf.size(); ++i)
    {
        num_all_meas += gnss_meas_buf[i].size();
        all_sat_states.push_back(sat_states(gnss_meas_buf[i], gnss_ephem_buf[i]));
    }
}

// 基于伪距进行粗略定位
bool GNSSVIInitializer::coarse_localization(Eigen::Matrix<double, 7, 1> &result)
{
    result.setZero();
    std::vector<ObsPtr> accum_obs; // 观测值
    std::vector<EphemBasePtr> accum_ephems; // 星历
    for (uint32_t i = 0; i < gnss_meas_buf.size(); ++i)
    {
        std::copy(gnss_meas_buf[i].begin(), gnss_meas_buf[i].end(), std::back_inserter(accum_obs));
        std::copy(gnss_ephem_buf[i].begin(), gnss_ephem_buf[i].end(), std::back_inserter(accum_ephems));
    }
    
    // 根据伪距进行定位，返回：ECEF 中的接收器位置和4个星座的四个时钟偏差 ；
    Eigen::Matrix<double, 7, 1> xyzt = psr_pos(accum_obs, accum_ephems, iono_params); 
    if (xyzt.topLeftCorner<3, 1>().norm() == 0)
    {
        std::cerr << "Failed to obtain a rough reference location.\n";
        return false;
    }

    for (uint32_t k = 0; k < 4; ++k)
    {
        if (fabs(xyzt(k+3)) < 1)
            xyzt(k+3) = 0;          // not observed yet
    }

    result = xyzt;
    return true;
}

// 偏航角对齐
// 这里的偏航角对齐用的是多普勒观测（因为相对伪距更准确）
// baseline：将【VIO的速度】与【GNSS的多普勒速度】的残差构建最小二乘，优化变量是 【偏航角】和【接收机钟差】
bool GNSSVIInitializer::yaw_alignment(const std::vector<Eigen::Vector3d> &local_vs, 
    const Eigen::Vector3d &rough_anchor_ecef, double &aligned_yaw, double &rcv_ddt)
{
    aligned_yaw = 0;
    rcv_ddt = 0;

    double est_yaw = 0;
    double est_rcv_ddt = 0;

    Eigen::Matrix3d rough_R_ecef_enu = ecef2rotation(rough_anchor_ecef); // ENU -> ECEF的旋转矩阵
    
    // 开始进行迭代求解
    uint32_t align_iter = 0;
    double align_dx_norm = 1.0;
    while (align_iter < MAX_ITERATION && align_dx_norm > CONVERGENCE_EPSILON)
    {
        Eigen::MatrixXd align_G(num_all_meas, 2);
        align_G.setZero();
        align_G.col(1).setOnes();
        Eigen::VectorXd align_b(num_all_meas);
        align_b.setZero();
        
        // local 到 ENU 的旋转矩阵；
        Eigen::Matrix3d align_R_enu_local(Eigen::AngleAxisd(est_yaw, Eigen::Vector3d::UnitZ()));
        Eigen::Matrix3d align_tmp_M;

        // 绕Z轴旋转的矩阵的导数，为后面矩阵构造使用，推导一下公式；
        align_tmp_M << -sin(est_yaw), -cos(est_yaw), 0,
                        cos(est_yaw), -sin(est_yaw), 0,
                        0       , 0        , 0;
        
        uint32_t align_counter = 0;
        for (uint32_t i = 0; i < gnss_meas_buf.size(); ++i)
        {
            // 局部VIO速度 转换到 ECEH坐标系下；
            Eigen::Matrix<double, 4, 1> ecef_vel_ddt;
            ecef_vel_ddt.head<3>() = rough_R_ecef_enu * align_R_enu_local * local_vs[i];
            ecef_vel_ddt(3) = est_rcv_ddt; // 接收机时间偏移；
            
            Eigen::VectorXd epoch_res;
            Eigen::MatrixXd epoch_J;
            
            // 多普勒残差计算；G，b计算；
            dopp_res(ecef_vel_ddt, rough_anchor_ecef, gnss_meas_buf[i], all_sat_states[i], epoch_res, epoch_J);
            align_b.segment(align_counter, gnss_meas_buf[i].size()) = epoch_res;
            align_G.block(align_counter, 0, gnss_meas_buf[i].size(), 1) = 
                epoch_J.leftCols(3)*rough_R_ecef_enu*align_tmp_M*local_vs[i];
            align_counter += gnss_meas_buf[i].size();
        }
        // 最小二乘求解
        Eigen::VectorXd dx = -(align_G.transpose()*align_G).inverse() * align_G.transpose() * align_b;
        est_yaw += dx(0);
        est_rcv_ddt += dx(1);
        align_dx_norm = dx.norm();
        ++ align_iter;
    }

    if (align_iter > MAX_ITERATION)
    {
        std::cerr << "Fail to initialize yaw offset.\n";
        return false;
    }

    // GNSS-VIO的偏航角赋值给全局变量aligned_yaw ；
    aligned_yaw = est_yaw;
    if (aligned_yaw > M_PI)
        aligned_yaw -= floor(est_yaw/(2.0*M_PI) + 0.5) * (2.0*M_PI);
    else if (aligned_yaw < -M_PI)
        aligned_yaw -=  ceil(est_yaw/(2.0*M_PI) - 0.5) * (2.0*M_PI);

    // 接收机时间偏差赋值全局变量rcv_ddt；
    rcv_ddt = est_rcv_ddt;

    return true;
}

bool GNSSVIInitializer::anchor_refinement(const std::vector<Eigen::Vector3d> &local_ps, 
    const double aligned_yaw, const double aligned_ddt, 
    const Eigen::Matrix<double, 7, 1> &rough_ecef_dt, Eigen::Matrix<double, 7, 1> &refined_ecef_dt)
{
    refined_ecef_dt.setZero();

    Eigen::Matrix3d aligned_R_enu_local(Eigen::AngleAxisd(aligned_yaw, Eigen::Vector3d::UnitZ()));

    // refine anchor point and receiver clock bias
    Eigen::Vector3d refine_anchor = rough_ecef_dt.head<3>();
    Eigen::Vector4d refine_dt = rough_ecef_dt.tail<4>();
    uint32_t refine_iter = 0;
    double refine_dx_norm = 1.0;
    std::vector<uint32_t> unobserved_sys;
    for (uint32_t k = 0; k < 4; ++k)
    {
        if (rough_ecef_dt(3+k) == 0)
            unobserved_sys.push_back(k);
    }

    while (refine_iter < MAX_ITERATION && refine_dx_norm > CONVERGENCE_EPSILON)
    {
        Eigen::MatrixXd refine_G(num_all_meas+unobserved_sys.size(), 7);
        Eigen::VectorXd refine_b(num_all_meas+unobserved_sys.size());
        refine_G.setZero();
        refine_b.setZero();
        uint32_t refine_counter = 0;
        Eigen::Matrix3d refine_R_ecef_enu = ecef2rotation(refine_anchor);
        Eigen::Matrix3d refine_R_ecef_local = refine_R_ecef_enu * aligned_R_enu_local;
        for (uint32_t i = 0; i < gnss_meas_buf.size(); ++i)
        {
            Eigen::Matrix<double, 7, 1> ecef_xyz_dt;
            ecef_xyz_dt.head<3>() = refine_R_ecef_local * local_ps[i] + refine_anchor;
            ecef_xyz_dt.tail<4>() = refine_dt + aligned_ddt * i * Eigen::Vector4d::Ones();

            Eigen::VectorXd epoch_res;
            Eigen::MatrixXd epoch_J;
            std::vector<Eigen::Vector2d> tmp_atmos_delay, tmp_sv_azel;
            psr_res(ecef_xyz_dt, gnss_meas_buf[i], all_sat_states[i], iono_params, 
                epoch_res, epoch_J, tmp_atmos_delay, tmp_sv_azel);
            refine_b.segment(refine_counter, gnss_meas_buf[i].size()) = epoch_res;
            refine_G.middleRows(refine_counter, gnss_meas_buf[i].size()) = epoch_J;
            refine_counter += gnss_meas_buf[i].size();
        }
        for (uint32_t k : unobserved_sys)
        {
            refine_b(refine_counter) = 0;
            refine_G(refine_counter, k+3) = 1.0;
            ++ refine_counter;
        }

        Eigen::VectorXd dx = -(refine_G.transpose()*refine_G).inverse() * refine_G.transpose() * refine_b;
        refine_anchor += dx.head<3>();
        refine_dt += dx.tail<4>();
        refine_dx_norm = dx.norm();
        ++ refine_iter;
    }

    if (refine_iter > MAX_ITERATION)
    {
        std::cerr << "Fail to perform anchor refinement.\n";
        return false;
    }

    refined_ecef_dt.head<3>() = refine_anchor;
    refined_ecef_dt.tail<4>() = refine_dt;

    return true;
}