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

std::shared_ptr<RDKit::ROMol> optimizeH(RDKit::ROMol& mol, double forceConst, bool addH) {
    // 创建MMFF力场
    std::shared_ptr<RDKit::ROMol> mol_tmp;
    if (addH) {
        mol_tmp.reset(RDKit::MolOps::addHs(mol, false, true));
    }
    else {
        mol_tmp = std::make_shared<RDKit::ROMol>(mol);
    }

    RDKit::MMFF::MMFFMolProperties mmffProps(*mol_tmp, "MMFF94s");
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
        RDKit::MMFF::constructForceField(*mol_tmp, &mmffProps));
    
    if (!ff) {
        throw std::runtime_error("Failed to construct force field");
    }
    
    // 对重原子添加位置约束（而不是完全固定）
    for (unsigned int i = 0; i < mol_tmp->getNumAtoms(); ++i) {
        if (mol_tmp->getAtomWithIdx(i)->getAtomicNum() > 1) {
            ForceFields::ContribPtr pc(
                new ForceFields::MMFF::PositionConstraintContrib(
                    ff.get(), i, 0.0, forceConst));  // 适中的约束强度
            ff->contribs().push_back(pc);
        }
    }
    
    // 最小化氢原子位置
    ff->minimize(1000);
    return mol_tmp;
}

std::shared_ptr<RDKit::ROMol> MiniFixAtomTor(RDKit::ROMol& mol, 
                    const std::vector<int>& fixatoms,  
                    const std::vector<int>& idx_,
                    double forceConst,
                    bool addH) {    
    std::shared_ptr<RDKit::ROMol> mol_tmp;
    if (addH) {
        mol_tmp.reset(RDKit::MolOps::addHs(mol, false, true));
    }
    else {
        mol_tmp = std::make_shared<RDKit::ROMol>(mol);
    }
    // 获取原子数量，并确保类型匹配
    unsigned int num_atoms = mol_tmp->getNumAtoms();
    
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
    RDKit::MMFF::MMFFMolProperties mmffProps(*mol_tmp, "MMFF94s");
    std::unique_ptr<ForceFields::ForceField> ff(
        RDKit::MMFF::constructForceField(*mol_tmp, &mmffProps, 1e9, -1, true));
    
    if (!ff) {
        throw std::runtime_error("Could not create MMFF94s force field");
    }

    mmffProps.setMMFFDielectricConstant(4.0);
    // 设置介电模型（1 表示常数介电模型）
    mmffProps.setMMFFDielectricModel(2);

    // 使用位置约束代替固定原子位置
    for (int fixid : fixatoms) {
        if (mol_tmp->getAtomWithIdx(fixid)->getAtomicNum() > 1) {
            // 添加位置约束而不是完全固定
            ForceFields::ContribPtr pc(
                new ForceFields::MMFF::PositionConstraintContrib(
                    ff.get(), fixid, 0.0, forceConst));  
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
                        
                        RDKit::Bond* bond1 = mol_tmp->getBondBetweenAtoms(idx_[i], idx_[j]);
                        RDKit::Bond* bond2 = mol_tmp->getBondBetweenAtoms(idx_[j], idx_[k]);
                        RDKit::Bond* bond3 = mol_tmp->getBondBetweenAtoms(idx_[k], idx_[l]);
                        
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
    
    ff->minimize(1000);
    return mol_tmp;
}