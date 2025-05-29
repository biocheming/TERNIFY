#include "align.hpp"
#include <iostream>
#include <unordered_map>
#include <memory>

// RDKit includes
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolTransforms/MolTransforms.h>


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


// 找出两个分子中对应的可旋转二面角
std::vector<std::pair<std::array<int, 4>, std::array<int, 4>>> 
findCorrespondingDihedrals(const RDKit::ROMol& mol1, 
                          const RDKit::ROMol& mol2, 
                          const RDKit::MatchVectType& match) {
    
    // 避免未使用参数警告
    (void)mol2;
    
    std::vector<std::pair<std::array<int, 4>, std::array<int, 4>>> dihedral_pairs;
    
    // 创建原子索引映射：mol1 -> mol2
    std::unordered_map<int, int> atom_map;
    for (const auto& pair : match) {
        atom_map[pair.second] = pair.first;  // mol1_idx -> mol2_idx
    }
    
    // 定义可旋转键的SMARTS模式
    std::string rotatable_smarts = "[!$(C(=[N,O,S])-!@[#7H,O,S])&!$([#7H,O,S]-!@C=[N,O,S])&!D1]-&!@[!D1]";
    std::unique_ptr<RDKit::ROMol> rot_pattern(RDKit::SmartsToMol(rotatable_smarts));
    
    if (!rot_pattern) {
        std::cerr << "Failed to parse rotatable bond SMARTS pattern" << std::endl;
        return dihedral_pairs;
    }
    
    // 在mol1中找到可旋转键
    std::vector<RDKit::MatchVectType> rot_matches1 = RDKit::SubstructMatch(mol1, *rot_pattern);
    
    for (const auto& rot_match : rot_matches1) {
        if (rot_match.size() >= 2) {
            int atom1_mol1 = rot_match[0].second;
            int atom2_mol1 = rot_match[1].second;
            
            // 检查这两个原子是否都在匹配中
            if (atom_map.find(atom1_mol1) != atom_map.end() && 
                atom_map.find(atom2_mol1) != atom_map.end()) {
                
                int atom1_mol2 = atom_map[atom1_mol1];
                int atom2_mol2 = atom_map[atom2_mol1];
                
                // 为mol1中的键找到二面角
                std::array<int, 4> dihedral_mol1 = {-1, -1, -1, -1};
                std::array<int, 4> dihedral_mol2 = {-1, -1, -1, -1};
                
                // 找到atom1的邻居（不是atom2）
                const RDKit::Atom* a1 = mol1.getAtomWithIdx(atom1_mol1);
                auto [nbrIdx1, endNbrs1] = mol1.getAtomNeighbors(a1);
                while (nbrIdx1 != endNbrs1) {
                    if (static_cast<int>(*nbrIdx1) != atom2_mol1) {
                        dihedral_mol1[0] = static_cast<int>(*nbrIdx1);
                        break;
                    }
                    ++nbrIdx1;
                }
                
                // 找到atom2的邻居（不是atom1）
                const RDKit::Atom* a2 = mol1.getAtomWithIdx(atom2_mol1);
                auto [nbrIdx2, endNbrs2] = mol1.getAtomNeighbors(a2);
                while (nbrIdx2 != endNbrs2) {
                    if (static_cast<int>(*nbrIdx2) != atom1_mol1) {
                        dihedral_mol1[3] = static_cast<int>(*nbrIdx2);
                        break;
                    }
                    ++nbrIdx2;
                }
                
                if (dihedral_mol1[0] != -1 && dihedral_mol1[3] != -1) {
                    dihedral_mol1[1] = atom1_mol1;
                    dihedral_mol1[2] = atom2_mol1;
                    
                    // 为mol2构建对应的二面角
                    if (atom_map.find(dihedral_mol1[0]) != atom_map.end() && 
                        atom_map.find(dihedral_mol1[3]) != atom_map.end()) {
                        
                        dihedral_mol2[0] = atom_map[dihedral_mol1[0]];
                        dihedral_mol2[1] = atom1_mol2;
                        dihedral_mol2[2] = atom2_mol2;
                        dihedral_mol2[3] = atom_map[dihedral_mol1[3]];
                        
                        dihedral_pairs.push_back({dihedral_mol1, dihedral_mol2});
                        
                        std::cout << "Found corresponding dihedral: "
                                  << "mol1[" << dihedral_mol1[0] << "-" << dihedral_mol1[1] 
                                  << "-" << dihedral_mol1[2] << "-" << dihedral_mol1[3] << "] <-> "
                                  << "mol2[" << dihedral_mol2[0] << "-" << dihedral_mol2[1] 
                                  << "-" << dihedral_mol2[2] << "-" << dihedral_mol2[3] << "]" << std::endl;
                    }
                }
            }
        }
    }
    
    return dihedral_pairs;
}

// 调整分子中的二面角使其与参考分子一致
void alignDihedrals(RDKit::ROMol& mol_to_adjust, 
                   const RDKit::ROMol& reference_mol,
                   const std::vector<std::pair<std::array<int, 4>, std::array<int, 4>>>& dihedral_pairs) {
    
    RDKit::Conformer& conf_adjust = mol_to_adjust.getConformer();
    const RDKit::Conformer& conf_ref = reference_mol.getConformer();
    
    for (const auto& pair : dihedral_pairs) {
        const auto& dihedral_adjust = pair.first;   // mol_to_adjust中的二面角
        const auto& dihedral_ref = pair.second;     // reference_mol中的二面角
        
        try {
            // 获取参考分子中的二面角值 - 使用 MolTransforms 而不是 RDKit::MolTransforms
            double ref_angle = MolTransforms::getDihedralDeg(conf_ref,
                                                           dihedral_ref[0],
                                                           dihedral_ref[1],
                                                           dihedral_ref[2],
                                                           dihedral_ref[3]);
            
            // 设置要调整分子中的二面角值 - 使用 MolTransforms 而不是 RDKit::MolTransforms
            MolTransforms::setDihedralDeg(conf_adjust,
                                        dihedral_adjust[0],
                                        dihedral_adjust[1],
                                        dihedral_adjust[2],
                                        dihedral_adjust[3],
                                        ref_angle);
            
            std::cout << "Aligned dihedral [" << dihedral_adjust[0] << "-" << dihedral_adjust[1] 
                      << "-" << dihedral_adjust[2] << "-" << dihedral_adjust[3] 
                      << "] to " << ref_angle << " degrees" << std::endl;
                      
        } catch (const std::exception& e) {
            std::cerr << "Error aligning dihedral: " << e.what() << std::endl;
        }
    }
}