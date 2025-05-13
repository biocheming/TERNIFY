#include "protein.hpp"
#include "parameters.hpp"
#include <fstream>
#include <string>
#include <algorithm>

Protein::Protein() = default;

void Protein::ReadProt(const std::string& filename, const std::array<double, 3>& translation) {
    std::ifstream file(filename);
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.substr(0, 4) == "ATOM") {
            // 检查是否为氢原子且是主链或A构象
            if (line[13] != 'H' && (line[16] == ' ' || line[16] == 'A')) {
                // 保存结构信息
                struct_data.push_back(line);
                
                // 解析坐标
                std::array<double, 3> coord;
                coord[0] = std::stod(line.substr(30, 8));
                coord[1] = std::stod(line.substr(38, 8));
                coord[2] = std::stod(line.substr(46, 8));
                coords.push_back(coord);
                
                // 获取残基和原子类型
                std::string res = line.substr(17, 3);
                if (res == "HIS") res = "HSD";
                
                std::string atom_type = line.substr(13, 3);
                // 去除atom_type前后的空格
                atom_type.erase(0, atom_type.find_first_not_of(" "));
                atom_type.erase(atom_type.find_last_not_of(" ") + 1);
                
                std::string atp = res + "_" + atom_type;
                
                // 检查原子类型并设置参数
                auto atom_it = atomtypes.find(atp);
                if (atom_it != atomtypes.end()) {
                    auto hbond_it = hbondtypes.find(atp);
                    if (hbond_it != hbondtypes.end()) {
                        if (hbondtypes[atp]== 2) {
                            para.push_back(std::make_tuple(std::nullopt, atomtypes[atp], 2));
                        } else { //hbondtypes[atp]== 3
                            para.push_back(std::make_tuple(std::nullopt, atomtypes[atp], 3));
                        }
                    } else {
                        para.push_back(std::make_tuple(std::nullopt, atomtypes[atp], 0));
                    }
                } else {
                    para.push_back(std::make_tuple(std::nullopt, 0.0, 0));
                }
            }
        }
        else if (line.substr(0, 3) == "END") {
            break;
        }
        else if (line.substr(0, 3) == "TER") {
            struct_data.push_back(line);
        }
    }
    
    // 应用平移
    for (auto& coord : coords) {
        for (int i = 0; i < 3; ++i) {
            coord[i] -= translation[i];
        }
    }
}