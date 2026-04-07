#include "pose_local_parameterization.h"

bool PoseLocalParameterization::Plus(const double *x, const double *delta, double *x_plus_delta) const
{
    Eigen::Map<const Eigen::Vector3d> _p(x); // 位置
    Eigen::Map<const Eigen::Quaterniond> _q(x + 3); // 旋转（四元数）

    Eigen::Map<const Eigen::Vector3d> dp(delta);  // 位置的增量

    Eigen::Quaterniond dq = Utility::deltaQ(Eigen::Map<const Eigen::Vector3d>(delta + 3));  // 旋转的增量

    Eigen::Map<Eigen::Vector3d> p(x_plus_delta); // 更新后的位置
    Eigen::Map<Eigen::Quaterniond> q(x_plus_delta + 3); // 更新后的旋转

    p = _p + dp; // 位置增量 
    q = (_q * dq).normalized(); // 旋转增量，使用四元数乘法，并归一化

    return true;
}
bool PoseLocalParameterization::ComputeJacobian(const double *x, double *jacobian) const
{
    Eigen::Map<Eigen::Matrix<double, 7, 6, Eigen::RowMajor>> j(jacobian);
    j.topRows<6>().setIdentity(); // 设置前 6 行为单位矩阵
    j.bottomRows<1>().setZero();  // 设置最后一行为零

    return true;
}