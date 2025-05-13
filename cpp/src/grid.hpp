

#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <tuple>
#include "common.hpp"
#include "parameters.hpp"
#include <optional>

typedef std::pair<GridValue, VolRegion> GridVR;
typedef std::vector<GridVR> Grids;
typedef std::pair<std::vector<std::vector<std::vector<GridValue>>>, VolRegion> GRID;
GRID Grid(const std::string& fpro, const VolRegion& pocket, int processes);

//一个格点GridValue和8个顶点VolRegion组成了GridVR
//若干个GridVR组成了Grids
//由Grids的索引和口袋VolRegion 组成了GRID