#include "minimize_h.hpp"
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <memory>
#include <vector>

std::shared_ptr<RDKit::ROMol> MinimizeH(const RDKit::ROMol& input_mol, double forceConst, bool addH_only) {
    
    if (addH_only) {
        // First create a copy of the input molecule without hydrogens
        std::shared_ptr<RDKit::ROMol> mol_no_h(RDKit::MolOps::removeHs(input_mol));
        // Then add hydrogens with 3D coordinates
        std::shared_ptr<RDKit::ROMol> mol(RDKit::MolOps::addHs(*mol_no_h, false, true));
        return mol;
    }
    
    // 添加氢原子并生成坐标
    std::shared_ptr<RDKit::ROMol> mol(RDKit::MolOps::addHs(input_mol, false, true));
    // 创建MMFF力场
    RDKit::MMFF::MMFFMolProperties mmffProps(*mol, "MMFF94s");
    if (!mmffProps.isValid()) {
        throw std::runtime_error("MMFF properties are invalid for this molecule.");
    }
    
    // 关闭VDW和静电相互作用
    mmffProps.setMMFFVdWTerm(false);
    mmffProps.setMMFFEleTerm(false);
    
    // 构建力场
    std::unique_ptr<ForceFields::ForceField> ff(
        RDKit::MMFF::constructForceField(*mol, &mmffProps));
    
    if (!ff) {
        throw std::runtime_error("Failed to construct force field");
    }
    
    // 对重原子添加位置约束（而不是完全固定）
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
        if (mol->getAtomWithIdx(i)->getAtomicNum() > 1) {
            ForceFields::ContribPtr pc(
                new ForceFields::MMFF::PositionConstraintContrib(
                    ff.get(), i, 0.0, forceConst));  // 适中的约束强度
            ff->contribs().push_back(pc);
        }
    }
    
    // 最小化氢原子位置
    ff->minimize(10000);
    
    return mol;
}

void optimizeH(RDKit::ROMol& mol, double forceConst) {
    // 创建MMFF力场
    RDKit::MMFF::MMFFMolProperties mmffProps(mol, "MMFF94s");
    if (!mmffProps.isValid()) {
        throw std::runtime_error("MMFF properties are invalid for this molecule.");
    }
    
    // 关闭VDW和静电相互作用
    //mmffProps.setMMFFVdWTerm(false);
    //mmffProps.setMMFFEleTerm(false);
    mmffProps.setMMFFDielectricConstant(4.0);

    // 设置介电模型（1 表示常数介电模型）
    mmffProps.setMMFFDielectricModel(2);
    
    // 构建力场
    std::unique_ptr<ForceFields::ForceField> ff(
        RDKit::MMFF::constructForceField(mol, &mmffProps));
    
    if (!ff) {
        throw std::runtime_error("Failed to construct force field");
    }
    
    // 对重原子添加位置约束（而不是完全固定）
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        if (mol.getAtomWithIdx(i)->getAtomicNum() > 1) {
            ForceFields::ContribPtr pc(
                new ForceFields::MMFF::PositionConstraintContrib(
                    ff.get(), i, 0.0, forceConst));  // 适中的约束强度
            ff->contribs().push_back(pc);
        }
    }
    
    // 最小化氢原子位置
    ff->minimize(1000);
}

void optimizeWithConstrAtoms(RDKit::ROMol& mol, const std::vector<int>& fixedAtoms) {
    // 初始化 MMFF 参数
    RDKit::MMFF::MMFFMolProperties mmffProps(mol, "MMFF94s");
    if (!mmffProps.isValid()) {
        std::cerr << "MMFF properties are invalid for this molecule." << std::endl;
        return;
    }

    std::unique_ptr<ForceFields::ForceField> ff(
        RDKit::MMFF::constructForceField(mol, &mmffProps, 1e9, -1, true));
    
    if (!ff) {
        std::cerr << "Failed to initialize MMFF force field." << std::endl;
        return;
    }

    // 对固定原子添加位置约束（而不是完全固定）
    for (size_t i = 0; i < fixedAtoms.size(); ++i) {
        int idx = fixedAtoms[i];
        if (idx >= 0 && static_cast<unsigned int>(idx) < mol.getNumAtoms() && 
            mol.getAtomWithIdx(idx)->getAtomicNum() > 1) {
            
            // 使用适中的力常数进行位置约束
            ForceFields::ContribPtr pc(
                new ForceFields::MMFF::PositionConstraintContrib(
                    ff.get(), idx, 0.0, 100.0));  // 适中的约束强度，允许小幅移动
            ff->contribs().push_back(pc);
        }
    }
    // 执行力场优化
    int maxIters = 500;  // 设置最大迭代次数
    ff->minimize(maxIters);
}


void optimizeWithFixedAtoms(RDKit::ROMol& mol, const std::vector<int>& fixedAtoms) {
    // 初始化 MMFF 参数
    RDKit::MMFF::MMFFMolProperties mmffProps(mol, "MMFF94s");
    if (!mmffProps.isValid()) {
        std::cerr << "MMFF properties are invalid for this molecule." << std::endl;
        return;
    }

    // 第一步：记录固定原子的初始坐标
    std::vector<RDGeom::Point3D> fixed_poses;
    RDKit::Conformer& conf = mol.getConformer();
    for (int idx : fixedAtoms) {
        if (idx >= 0 && static_cast<unsigned int>(idx) < mol.getNumAtoms()) {
            fixed_poses.push_back(conf.getAtomPos(idx));
        }
    }

    // 第一步：使用位置约束进行温和优化
    std::unique_ptr<ForceFields::ForceField> ff1(
        RDKit::MMFF::constructForceField(mol, &mmffProps, 1e9, -1, true));
    
    if (!ff1) {
        std::cerr << "Failed to initialize MMFF force field for constrained optimization." << std::endl;
        return;
    }

    // 对固定原子添加位置约束（而不是完全固定）
    for (size_t i = 0; i < fixedAtoms.size(); ++i) {
        int idx = fixedAtoms[i];
        if (idx >= 0 && static_cast<unsigned int>(idx) < mol.getNumAtoms() && 
            mol.getAtomWithIdx(idx)->getAtomicNum() > 1) {
            
            // 使用适中的力常数进行位置约束
            ForceFields::ContribPtr pc(
                new ForceFields::MMFF::PositionConstraintContrib(
                    ff1.get(), idx, 0.0, 100.0));  // 适中的约束强度，允许小幅移动
            ff1->contribs().push_back(pc);
        }
    }

    // 执行第一轮优化（位置约束优化）
    int maxIters1 = 500;
    ff1->minimize(maxIters1);

    // 第二步：恢复固定原子的原始坐标
    for (size_t i = 0; i < fixedAtoms.size(); ++i) {
        int idx = fixedAtoms[i];
        if (idx >= 0 && static_cast<unsigned int>(idx) < mol.getNumAtoms()) {
            conf.setAtomPos(idx, fixed_poses[i]);
        }
    }

    // 第二步：完全固定这些原子进行最终优化
    std::unique_ptr<ForceFields::ForceField> ff2(
        RDKit::MMFF::constructForceField(mol, &mmffProps, 1e9, -1, true));
    
    if (!ff2) {
        std::cerr << "Failed to initialize MMFF force field for fixed optimization." << std::endl;
        return;
    }

    // 将指定的原子完全固定
    for (int idx : fixedAtoms) {
        if (idx >= 0 && static_cast<unsigned int>(idx) < mol.getNumAtoms() && 
            mol.getAtomWithIdx(idx)->getAtomicNum() > 1) {
            ff2->fixedPoints().push_back(idx);
        }
    }

    // 执行第二轮优化（固定原子优化）
    int maxIters2 = 500;
    ff2->minimize(maxIters2);
}

void MiniFixAtomTor(RDKit::ROMol& mol, 
                    const std::vector<int>& fixatoms,  
                    const std::vector<int>& idx_) {    
    // 获取原子数量，并确保类型匹配
    unsigned int num_atoms = mol.getNumAtoms();
    
    // 检查输入的原子索引是否有效
    for (int idx : fixatoms) {
        if (idx < 0 || static_cast<unsigned int>(idx) >= num_atoms) {
            throw std::runtime_error("Invalid fixed atom index: " + std::to_string(idx));
        }
    }
    for (int idx : idx_) {
        if (idx < 0 || static_cast<unsigned int>(idx) >= num_atoms) {
            throw std::runtime_error("Invalid constrained atom index: " + std::to_string(idx));
        }
    }

    // 创建力场
    RDKit::MMFF::MMFFMolProperties mmffProps(mol, "MMFF94s");
    std::unique_ptr<ForceFields::ForceField> ff(
        RDKit::MMFF::constructForceField(mol, &mmffProps, 1e9, -1, true));
    
    if (!ff) {
        throw std::runtime_error("Could not create MMFF94s force field");
    }
    // 关闭VDW和静电相互作用
    mmffProps.setMMFFVdWTerm(false);
    mmffProps.setMMFFEleTerm(false);
/*    
    // 固定原子位置
    for (int fixid : fixatoms) {
        if (mol.getAtomWithIdx(fixid)->getAtomicNum() > 1) {
            ff->fixedPoints().push_back(fixid);
        }
    }
*/
    // 使用位置约束代替固定原子位置
    for (int fixid : fixatoms) {
        if (mol.getAtomWithIdx(fixid)->getAtomicNum() > 1) {
            // 添加位置约束而不是完全固定
            ForceFields::ContribPtr pc(
                new ForceFields::MMFF::PositionConstraintContrib(
                    ff.get(), fixid, 0.0, 100.0));  // 使用20.0的力常数进行位置约束
            ff->contribs().push_back(pc);
        }
    }

    // 添加二面角约束
    for (size_t i = 0; i < idx_.size(); ++i) {
        for (size_t j = 0; j < idx_.size(); ++j) {
            for (size_t k = 0; k < idx_.size(); ++k) {
                for (size_t l = 0; l < idx_.size(); ++l) {
                    if (i != j && j != k && k != l && 
                        i != k && j != l && i != l) {
                        
                        RDKit::Bond* bond1 = mol.getBondBetweenAtoms(idx_[i], idx_[j]);
                        RDKit::Bond* bond2 = mol.getBondBetweenAtoms(idx_[j], idx_[k]);
                        RDKit::Bond* bond3 = mol.getBondBetweenAtoms(idx_[k], idx_[l]);
                        
                        if (bond1 && bond2 && bond3) {
                            auto *torsionConstraint = new ForceFields::MMFF::TorsionConstraintContrib(
                                ff.get(),
                                idx_[i], 
                                idx_[j], 
                                idx_[k], 
                                idx_[l],
                                true, //relative dihedral angle change range: [-0.5, 0.5]
                                -0.5,
                                0.5,
                                0.1 // force constant, never greater than 0.5 kcal/mol
                            );
                            ff->contribs().push_back(ForceFields::ContribPtr(torsionConstraint));
                        }
                    }
                }
            }
        }
    }
    
    ff->minimize(100);
}