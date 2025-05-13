#ifndef LIST_LOW_DIHE_HPP
#define LIST_LOW_DIHE_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm> // 用于 std::sort

class DihedralEnergy {
private:
    double v1, v2, v3; // 力场参数

public:
    // 构造函数
    DihedralEnergy(double v1, double v2, double v3) 
        : v1(v1), v2(v2), v3(v3) {}

    // 势能函数
    double potentialEnergy(double x) const {
        return 0.5 * (v1 * (1.0 + std::cos(x)) +
                      v2 * (1.0 - std::cos(2.0 * x)) +
                      v3 * (1.0 + std::cos(3.0 * x)));
    }

    // 势能函数的一阶导数
    double firstDerivative(double x) const {
        return -0.5 * (v1 * std::sin(x) - 2 * v2 * std::sin(2 * x) + 3 * v3 * std::sin(3 * x));
    }

    // 势能函数的二阶导数
    double secondDerivative(double x) const {
        return -0.5 * (v1 * std::cos(x) - 4 * v2 * std::cos(2.0 * x) + 9 * v3 * std::cos(3.0 * x));
    }

    // 寻找极小值并返回结果
    std::vector<double> findLocalMinima(double start, double end, double step) const {
        std::vector<double> minima; // 存储极小值点
        for (double x = start; x <= end; x += step) {
            double f_prime = firstDerivative(x);
            double f_double_prime = secondDerivative(x);
            if (std::abs(f_prime) < 1e-3 && f_double_prime > 0) { // 判断一阶导数是否接近0
                bool is_duplicate = false;
                for (double m : minima) {
                    if (std::abs(m - x) < 0.1) { // 去重
                        is_duplicate = true;
                        break;
                    }
                }
                if (!is_duplicate) {
                    minima.push_back(x);
                }
            }
        }
        // 如果没有找到极小值点，至少返回一个点（能量最低的点）
        if (minima.empty())
            minima.push_back(std::numeric_limits<double>::min());
        return minima;
    }
    //均匀插值到36个点
    std::vector<double> augmentMinima(const std::vector<double>& minima, double start, double end, size_t targetSize) const {
        std::vector<double> augmentedMinima = minima;
        std::vector<std::pair<double, double>> all_points;

        // 首先添加所有极小值点
        for (double x : augmentedMinima) {
            all_points.emplace_back(x, potentialEnergy(x));
        }

        // 在每个极小值点附近采样
        size_t current_points = all_points.size();
        size_t points_needed = targetSize - current_points;
        size_t points_per_minimum = points_needed / (current_points > 0 ? current_points : 1) + 1;

        for (double x : augmentedMinima) {
            // 在每个极小值点两侧逐渐扩展采样
            for (size_t i = 1; i <= points_per_minimum; ++i) {
                // 使用二次函数增加采样间距，避免点过于集中
                double offset = i * i * 0.025;
                
                // 左侧点
                double new_x1 = x - offset;
                if (new_x1 >= start) {
                    bool too_close = false;
                    for (const auto& point : all_points) {
                        if (std::abs(new_x1 - point.first) < 0.1) {
                            too_close = true;
                            break;
                        }
                    }
                    if (!too_close) {
                        all_points.emplace_back(new_x1, potentialEnergy(new_x1));
                    }
                }

                // 右侧点
                double new_x2 = x + offset;
                if (new_x2 <= end) {
                    bool too_close = false;
                    for (const auto& point : all_points) {
                        if (std::abs(new_x2 - point.first) < 0.1) {
                            too_close = true;
                            break;
                        }
                    }
                    if (!too_close) {
                        all_points.emplace_back(new_x2, potentialEnergy(new_x2));
                    }
                }
            }
        }

        // 如果还是不够点，在整个范围内均匀添加
        if (all_points.size() < targetSize) {
            double step = (end - start) / (targetSize * 2);
            for (double x = start; x <= end && all_points.size() < targetSize * 2; x += step) {
                bool too_close = false;
                for (const auto& point : all_points) {
                    if (std::abs(x - point.first) < 0.1) {
                        too_close = true;
                        break;
                    }
                }
                if (!too_close) {
                    all_points.emplace_back(x, potentialEnergy(x));
                }
            }
        }

        // 按能量排序
        std::sort(all_points.begin(), all_points.end(),
                 [](const auto& a, const auto& b) { return a.second < b.second; });

        // 选择能量最低的targetSize个点
        std::vector<double> result;
        result.reserve(targetSize);
        
        for (const auto& point : all_points) {
            if (result.size() >= targetSize) break;
            result.push_back(point.first * 180.0 / M_PI);
        }

        // 确保有足够的点
        while (result.size() < targetSize) {
            // 在最大间隙处添加点
            double max_gap = 0.0;
            double insert_val = 0.0;
            for (size_t i = 0; i < result.size() - 1; ++i) {
                double gap = result[i + 1] - result[i];
                if (gap > max_gap) {
                    max_gap = gap;
                    insert_val = (result[i] + result[i + 1]) / 2.0;
                }
            }
            result.push_back(insert_val);
            std::sort(result.begin(), result.end());
        }

        assert(result.size() == targetSize);
        return result;
    }

};

#endif // LIST_LOW_DIHE_HPP
