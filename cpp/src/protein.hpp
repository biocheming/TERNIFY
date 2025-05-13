#include <vector>
#include <array>
#include <string>
#include <tuple>
#include <optional>

class Protein {
public:
    Protein();
    
    // 读取蛋白质文件
    void ReadProt(const std::string& filename, const std::array<double, 3>& translation);

    // 成员变量
    std::vector<std::array<double, 3>> coords;
    std::vector<std::string> struct_data;
    std::vector< std::tuple<std::optional<int> , double, std::optional<int> > > para;  
};