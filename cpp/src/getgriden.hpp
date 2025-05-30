#pragma once
#include <vector>
#include <array>
#include <tuple>
#include "grid.hpp"

// 计算网格能量
double GetGridEn(
    const GRID& grid,
    const Coords& coords,
    const std::vector<IQHb>& q);
