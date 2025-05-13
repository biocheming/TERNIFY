#include "getgriden.hpp"
#include "parameters.hpp"

double GetGridEn(
    const GRID& grid,
    const Coords& coords,
    const std::vector<IQHb>& q) {
    
    double vdw = 0.0;
    double ele = 0.0;
    double hb = 0.0;

    const auto& grid_values = grid.first;
    const auto& grid_bounds = grid.second;
    
    for (size_t i = 0; i < coords.size(); ++i) {
        const auto& coor = coords[i];
        
        // 检查坐标是否在网格范围内
        if (coor[0] >= grid_bounds[0][0] && coor[0] <= grid_bounds[0][1] &&
            coor[1] >= grid_bounds[1][0] && coor[1] <= grid_bounds[1][1] &&
            coor[2] >= grid_bounds[2][0] && coor[2] <= grid_bounds[2][1]) {
            
            // 计算网格索引
            std::array<int, 3> loc;
            for (int j = 0; j < 3; ++j) {
                loc[j] = static_cast<int>((coor[j] - grid_bounds[j][0]) / paras["grid_space"]);
            }
            
            // 获取体素值
            const auto& voxel = grid_values[loc[0]][loc[1]][loc[2]];
            
            // 计算各项能量贡献
            vdw += voxel.vdw_or_clash;
            ele += voxel.elec * std::get<1>(q[i]);
            
            // 氢键计算
            auto hd_type = std::get<2>(q[i]);
            if (hd_type.has_value()) {
                if (hd_type.value() == 2 && voxel.hd_donor.has_value()) { //体素记录的氢键角色要与移动中的角色保持一致
                    hb += voxel.hd_acceptor.value();
                } else if (hd_type.value() == 3 && voxel.hd_acceptor.has_value()) {
                    hb += voxel.hd_donor.value();
                }
            }
        }
    }
    
    return vdw + ele + hb;
}

/*
GridEnergyResult GetGridEn(
    const GRID& grid,
    const Coords& coords,
    const std::vector<IQHb>& q,
    bool calc_derivs = false) {
    
    GridEnergyResult result{0.0};
    if (calc_derivs) {
        result.derivatives.resize(coords.size(), {0.0, 0.0, 0.0});
    }

    const auto& grid_values = grid.first;
    const auto& grid_bounds = grid.second;
    
    for (size_t i = 0; i < coords.size(); ++i) {
        const auto& coor = coords[i];
        
        // 检查坐标是否在网格范围内
        if (coor[0] >= grid_bounds[0][0] && coor[0] <= grid_bounds[0][1] &&
            coor[1] >= grid_bounds[1][0] && coor[1] <= grid_bounds[1][1] &&
            coor[2] >= grid_bounds[2][0] && coor[2] <= grid_bounds[2][1]) {
            
            // 计算网格索引
            std::array<int, 3> loc;
            for (int j = 0; j < 3; ++j) {
                loc[j] = static_cast<int>((coor[j] - grid_bounds[j][0]) / paras["grid_space"]);
            }
            
            // 获取体素值
            const auto& voxel = grid_values[loc[0]][loc[1]][loc[2]];
            
            // 计算总能量
            double total_energy = voxel.vdw_or_clash;
            total_energy += voxel.elec * std::get<1>(q[i]);
            
            // 氢键计算
            auto hd_type = std::get<2>(q[i]);
            if (hd_type.has_value()) {
                if (hd_type.value() == 2 && voxel.hd_donor.has_value()) {
                    total_energy += voxel.hd_acceptor.value();
                } else if (hd_type.value() == 3 && voxel.hd_acceptor.has_value()) {
                    total_energy += voxel.hd_donor.value();
                }
            }
            
            result.energy += total_energy;

            // 如果需要计算导数
            if (calc_derivs) {
                // 对每个维度计算有限差分导数
                for (int j = 0; j < 3; ++j) {
                    if (loc[j] > 0 && loc[j] < static_cast<int>(grid_values.size())-1) {
                        // 获取相邻格点的能量
                        std::array<int, 3> loc_plus = loc;
                        std::array<int, 3> loc_minus = loc;
                        loc_plus[j]++;
                        loc_minus[j]--;

                        const auto& voxel_plus = grid_values[loc_plus[0]][loc_plus[1]][loc_plus[2]];
                        const auto& voxel_minus = grid_values[loc_minus[0]][loc_minus[1]][loc_minus[2]];

                        // 计算相邻格点的总能量
                        double e_plus = voxel_plus.vdw_or_clash + voxel_plus.elec * std::get<1>(q[i]);
                        double e_minus = voxel_minus.vdw_or_clash + voxel_minus.elec * std::get<1>(q[i]);

                        // 添加氢键能量
                        if (hd_type.has_value()) {
                            if (hd_type.value() == 2) {
                                if (voxel_plus.hd_donor.has_value()) 
                                    e_plus += voxel_plus.hd_acceptor.value();
                                if (voxel_minus.hd_donor.has_value()) 
                                    e_minus += voxel_minus.hd_acceptor.value();
                            } else if (hd_type.value() == 3) {
                                if (voxel_plus.hd_acceptor.has_value()) 
                                    e_plus += voxel_plus.hd_donor.value();
                                if (voxel_minus.hd_acceptor.has_value()) 
                                    e_minus += voxel_minus.hd_donor.value();
                            }
                        }

                        // 计算中心差分导数
                        result.derivatives[i][j] = (e_plus - e_minus) / (2.0 * paras["grid_space"]);
                    }
                }
            }
        }
    }
    
    return result;
}

*/