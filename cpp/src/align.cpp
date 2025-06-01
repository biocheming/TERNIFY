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
    
    // Python: covar = np.dot(coord_var.T, coord_ref)
    // This means coord_var.transpose() * coord_ref
    MatrixXd coord_var_matrix(coord_var.size(), 3);
    MatrixXd coord_ref_matrix(coord_ref.size(), 3);
    
    for (size_t i = 0; i < coord_var.size(); ++i) {
        coord_var_matrix.row(i) = Map<const Vector3d>(coord_var[i].data()).transpose();
        coord_ref_matrix.row(i) = Map<const Vector3d>(coord_ref[i].data()).transpose();
    }
    
    // Python: covar = np.dot(coord_var.T, coord_ref)
    MatrixXd covar = coord_var_matrix.transpose() * coord_ref_matrix;

    // SVD分解
    JacobiSVD<MatrixXd> svd(covar, ComputeFullU | ComputeFullV);
    Matrix3d u = svd.matrixU();
    Matrix3d vh = svd.matrixV().transpose();  // Note: V^T, not V

    // Python: d = (np.linalg.det(u) * np.linalg.det(vh)) < 0.0
    // Python: if d: u[:, -1] = -u[:, -1]
    double det = u.determinant() * vh.determinant();
    if (det < 0.0) {
        u.col(2) = -u.col(2);  // Flip the last column of u
    }
    
    // Python: R = np.dot(u, vh)
    return u * vh;
}

std::vector<std::array<double, 3>> Align(const std::vector<std::array<double, 3>>& coords,
                                         const std::vector<std::array<double, 3>>& coord_subs_var,
                                         const std::vector<std::array<double, 3>>& coord_ref) {
    // Python: center = coord_ref.mean(axis=0)
    Vector3d center = Vector3d::Zero();
    for (const auto& point : coord_ref) {
        center += Map<const Vector3d>(point.data());
    }
    center /= coord_ref.size();

    // Python: coord_ref = coord_ref - center
    std::vector<std::array<double, 3>> centered_ref;
    for (const auto& point : coord_ref) {
        Vector3d centered = Map<const Vector3d>(point.data()) - center;
        centered_ref.push_back({centered.x(), centered.y(), centered.z()});
    }

    // Python: R = Kabsch(coord_subs_var, coord_ref)
    // Note: Python version does NOT center coord_subs_var in Align function
    Matrix3d R = Kabsch(coord_subs_var, centered_ref);

    // Python: coords = np.dot(coords, R) + center
    // This means coords * R, not R * coords!
    std::vector<std::array<double, 3>> aligned_coords;
    for (const auto& point : coords) {
        Vector3d coord_point = Map<const Vector3d>(point.data());
        // Python's np.dot(coords, R) when coords is a row vector means coords * R
        // Which is equivalent to R.transpose() * coords when treating coords as column vector
        Vector3d rotated = R.transpose() * coord_point;
        Vector3d transformed = rotated + center;
        aligned_coords.push_back({transformed.x(), transformed.y(), transformed.z()});
    }

    return aligned_coords;
}

std::vector<std::array<double, 3>> Align2(const std::vector<std::array<double, 3>>& coords,
                                          const std::vector<std::array<double, 3>>& coord_ref,
                                          const std::vector<std::array<double, 3>>& coord_subs_var,
                                          const std::array<double, 3>& translation) {
    // Python: center = coord_subs_var.mean(axis=0)
    Vector3d center = Vector3d::Zero();
    for (const auto& point : coord_subs_var) {
        center += Map<const Vector3d>(point.data());
    }
    center /= coord_subs_var.size();

    // Python: coord_subs_var = coord_subs_var - center
    std::vector<std::array<double, 3>> centered_subs_var;
    for (const auto& point : coord_subs_var) {
        Vector3d centered = Map<const Vector3d>(point.data()) - center;
        centered_subs_var.push_back({centered.x(), centered.y(), centered.z()});
    }

    // Python: R = Kabsch(coord_subs_var, coord_ref)
    Matrix3d R = Kabsch(centered_subs_var, coord_ref);
    Vector3d trans_vec = Map<const Vector3d>(translation.data());
    
    // Python: coords = coords - center; coords = np.dot(coords, R) + translation
    std::vector<std::array<double, 3>> aligned_coords;
    for (const auto& point : coords) {
        // Python: coords = coords - center
        Vector3d centered = Map<const Vector3d>(point.data()) - center;
        // Python: coords = np.dot(coords, R) + translation
        // Python's np.dot(coords, R) means coords * R, which is R.transpose() * coords
        Vector3d rotated = R.transpose() * centered;
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