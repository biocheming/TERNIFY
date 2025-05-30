#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <vector>
#include <array>
#include <unordered_map>
#include <memory>
#include <Eigen/Dense>

// RDKit includes
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/Substruct/SubstructMatch.h>


using namespace Eigen;

// Function declarations

// Compute the Kabsch rotation matrix
Matrix3d Kabsch(const std::vector<std::array<double, 3>>& coord_var, const std::vector<std::array<double, 3>>& coord_ref);

// Align coordinates
std::vector<std::array<double, 3>> Align(const std::vector<std::array<double, 3>>& coords,
                                          const std::vector<std::array<double, 3>>& coord_subs_var,
                                          const std::vector<std::array<double, 3>>& coord_ref);

// Align coordinates with a specific translation
std::vector<std::array<double, 3>> Align2(const std::vector<std::array<double, 3>>& coords,
                                            const std::vector<std::array<double, 3>>& coord_ref,
                                            const std::vector<std::array<double, 3>>& coord_subs_var,
                                            const std::array<double, 3>& translation);

// 找出两个分子中对应的可旋转二面角
std::vector<std::pair<std::array<int, 4>, std::array<int, 4>>> 
findCorrespondingDihedrals(const RDKit::ROMol& mol1, 
                          const RDKit::ROMol& mol2, 
                          const RDKit::MatchVectType& match);

// 调整分子中的二面角使其与参考分子一致
void alignDihedrals(RDKit::ROMol& mol_to_adjust, 
                   const RDKit::ROMol& reference_mol,
                   const std::vector<std::pair<std::array<int, 4>, std::array<int, 4>>>& dihedral_pairs);

#endif // ALIGN_HPP