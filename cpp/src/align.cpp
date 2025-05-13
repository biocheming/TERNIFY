#include "align.hpp"
#include <iostream>
// Function definitions

Matrix3d Kabsch(const std::vector<std::array<double, 3>>& coord_var, 
                const std::vector<std::array<double, 3>>& coord_ref) {
    // check if the number of atoms is equal
    if (coord_var.size() != coord_ref.size()||
        coord_var.empty() || coord_ref.empty()) {
        throw std::invalid_argument("Kabsch Search Error: the number of atoms is not equal or empty.");
        exit(EXIT_FAILURE);
    }
    // 计算协方差矩阵
    MatrixXd covar = MatrixXd::Zero(3, 3);
    for (size_t i = 0; i < coord_var.size(); ++i) {
        covar += Map<const Vector3d>(coord_var[i].data()) * 
                 Map<const Vector3d>(coord_ref[i].data()).transpose();
    }

    // SVD分解
    JacobiSVD<MatrixXd> svd(covar, ComputeFullU | ComputeFullV);
    Matrix3d u = svd.matrixU();
    Matrix3d v = svd.matrixV();

    // 检查是否需要处理反射
    Matrix3d i = Matrix3d::Identity();
    double det = (u * v.transpose()).determinant();
    i(2,2) = std::copysign(1.0, det);
    // 计算最优旋转矩阵
    return u * i * v.transpose();
}

std::vector<std::array<double, 3>> Align(const std::vector<std::array<double, 3>>& coords,
                                         const std::vector<std::array<double, 3>>& coord_subs_var,
                                         const std::vector<std::array<double, 3>>& coord_ref) {
    // 计算质心
    Vector3d centroid_subs_var = Vector3d::Zero();
    Vector3d centroid_ref = Vector3d::Zero();

    for (const auto& point : coord_subs_var) {
        centroid_subs_var += Map<const Vector3d>(point.data());
    }
    centroid_subs_var /= coord_subs_var.size();

    for (const auto& point : coord_ref) {
        centroid_ref += Map<const Vector3d>(point.data());
    }
    centroid_ref /= coord_ref.size();

    // 将coord_subs_var和coord_ref移到原点
    std::vector<std::array<double, 3>> centered_subs_var, centered_ref;
    for (const auto& point : coord_subs_var) {
        Vector3d centered = Map<const Vector3d>(point.data()) - centroid_subs_var;
        centered_subs_var.push_back({centered.x(), centered.y(), centered.z()});
    }
    for (const auto& point : coord_ref) {
        Vector3d centered = Map<const Vector3d>(point.data()) - centroid_ref;
        centered_ref.push_back({centered.x(), centered.y(), centered.z()});
    }

    // 计算旋转矩阵（从coord_subs_var到coord_ref）
    Matrix3d R = Kabsch(centered_subs_var, centered_ref);

    // 计算平移向量
    //Vector3d translation = centroid_ref - centroid_subs_var;

    // 将蛋白质坐标相对于coord_subs_var进行变换
    std::vector<std::array<double, 3>> aligned_coords;
    for (const auto& point : coords) {
        // 使用旋转矩阵的转置（逆）来保持相对方向
        Vector3d transformed = R.transpose() * (Map<const Vector3d>(point.data()) - centroid_subs_var) + centroid_ref;
        aligned_coords.push_back({transformed.x(), transformed.y(), transformed.z()});
    }

    return aligned_coords;
}

std::vector<std::array<double, 3>> Align2(const std::vector<std::array<double, 3>>& coords,
                                          const std::vector<std::array<double, 3>>& coord_subs_var,
                                          const std::vector<std::array<double, 3>>& coord_ref,
                                          const std::array<double, 3>& translation) {
    // 计算原点位置的质心
    Vector3d center = Vector3d::Zero();
    for (const auto& point : coord_subs_var) {
        center += Map<const Vector3d>(point.data());
    }
    center /= coord_subs_var.size();

    // 将coord_subs_var相对于质心
    std::vector<std::array<double, 3>> centered_subs_var;
    for (const auto& point : coord_subs_var) {
        Vector3d centered = Map<const Vector3d>(point.data()) - center;
        centered_subs_var.push_back({centered.x(), centered.y(), centered.z()});
    }

    // 计算从原点位置到当前位置的旋转矩阵
    Matrix3d R = Kabsch(centered_subs_var, coord_ref);
    Vector3d trans_vec = Map<const Vector3d>(translation.data());
    
    // 对coords应用相同的变换
    std::vector<std::array<double, 3>> aligned_coords;
    for (const auto& point : coords) {
        // 1. 相对于原点位置的质心
        Vector3d centered = Map<const Vector3d>(point.data()) - center;
        // 2. 应用旋转（左乘R）
        Vector3d rotated = R * centered;  // 修改：改为左乘
        // 3. 加上translation
        Vector3d transformed = rotated + trans_vec;
        aligned_coords.push_back({transformed.x(), transformed.y(), transformed.z()});
    }

    return aligned_coords;
}
