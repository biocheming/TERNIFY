#pragma once
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>  
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/DistanceConstraint.h>
#include <ForceField/MMFF/PositionConstraint.h>
#include <ForceField/MMFF/AngleConstraint.h>
#include <ForceField/MMFF/TorsionConstraint.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <memory>

// 返回智能指针以管理内存
std::shared_ptr<RDKit::ROMol> MinimizeH(const RDKit::ROMol& input_mol, double forceConst = 1.0e9, bool addH_only = false);
void optimizeH(RDKit::ROMol& mol, double forceConst = 1.0e9); 
void optimizeWithFixedAtoms(RDKit::ROMol& mol, const std::vector<int>& fixedAtoms);
void optimizeWithConstrAtoms(RDKit::ROMol& mol, const std::vector<int>& fixedAtoms);
void MiniFixAtomTor(RDKit::ROMol& mol, 
                    const std::vector<int>& fixatoms,  // 固定位置的原子
                    const std::vector<int>& idx_);    // 内部几何构型需要固定的原子组   