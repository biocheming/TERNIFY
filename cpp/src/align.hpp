#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <vector>
#include <array>
#include <Eigen/Dense>

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
                                            const std::vector<std::array<double, 3>>& coord_subs_var,
                                            const std::vector<std::array<double, 3>>& coord_ref,
                                            const std::array<double, 3>& translation);

#endif // ALIGN_HPP