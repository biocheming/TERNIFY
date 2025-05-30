#include "getgriden.hpp"
#include "parameters.hpp"
#include <iostream>

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
            double atom_vdw = voxel.vdw_or_clash;
            double atom_ele = voxel.elec * std::get<1>(q[i]);
            double atom_hb = 0.0;
            
            // 氢键计算
            auto hd_type = std::get<2>(q[i]);
            if (hd_type.has_value()) {
                if (hd_type.value() == 2 && voxel.hd_donor.has_value()) { // acceptor原子访问donor网格能量
                    atom_hb = voxel.hd_donor.value();
                } else if (hd_type.value() == 3 && voxel.hd_acceptor.has_value()) { // donor原子访问acceptor网格能量
                    atom_hb = voxel.hd_acceptor.value();
                }
            }
            
            vdw += atom_vdw;
            ele += atom_ele;
            hb += atom_hb;
        }
    }

    return vdw + ele + hb;
}

