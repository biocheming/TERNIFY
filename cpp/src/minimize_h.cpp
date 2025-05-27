#include "minimize_h.hpp"
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <memory>
#include <vector>

std::shared_ptr<RDKit::ROMol> MinimizeH(const RDKit::ROMol& input_mol) {
    
    // 添加氢原子并生成坐标
    std::shared_ptr<RDKit::ROMol> mol(RDKit::MolOps::addHs(input_mol, false, true));
    
    // 创建MMFF力场
    RDKit::MMFF::MMFFMolProperties mmffProps(*mol);
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
                    ff.get(), i, 0.0, 1.0e5));  // 适中的约束强度
            ff->contribs().push_back(pc);
        }
    }
    
    // 最小化氢原子位置
    ff->minimize(10000);
    
    return mol;
}

void optimizeWithFixedAtoms(RDKit::ROMol& mol, const std::vector<int>& fixedAtoms) {
    // 初始化 MMFF 参数
    RDKit::MMFF::MMFFMolProperties mmffProps(mol);
    if (!mmffProps.isValid()) {
        std::cerr << "MMFF properties are invalid for this molecule." << std::endl;
        return;
    }

    std::unique_ptr<ForceFields::ForceField> ff(
        RDKit::MMFF::constructForceField(mol, &mmffProps, 1e6, -1, true));
    
    if (!ff) {
        std::cerr << "Failed to initialize MMFF force field." << std::endl;
        return;
    }

    // 将指定的原子固定
    for (int idx : fixedAtoms) {
        if (mol.getAtomWithIdx(idx)->getAtomicNum() > 1) {
            ff->fixedPoints().push_back(idx);
        }
    }
    // 执行力场优化
    int maxIters = 10000;  // 设置最大迭代次数
    ff->minimize(maxIters);
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
    RDKit::MMFF::MMFFMolProperties mmffProps(mol);
    std::unique_ptr<ForceFields::ForceField> ff(
        RDKit::MMFF::constructForceField(mol, &mmffProps, 1e6, -1, true));
    
    if (!ff) {
        throw std::runtime_error("Could not create MMFF94s force field");
    }
    
    // 固定原子位置
    for (int fixid : fixatoms) {
        if (mol.getAtomWithIdx(fixid)->getAtomicNum() > 1) {
            ff->fixedPoints().push_back(fixid);
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
    
    ff->minimize(5000);
}