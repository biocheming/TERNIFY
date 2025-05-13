#include "grid.hpp"
#include <fstream>
#include <cmath>
#include <thread>
#include <future>
#include <vector>
#include <mutex>

#include <iostream>


std::array<double, 2> CalcEn(const std::vector<TerAtom>& atoms, 
                              const std::array<int, 3>& coor,
                              const std::array<std::array<double, 2>, 3>& pocket,
                              GridValue& value) {
    if (value.vdw_or_clash == 0 ) { // Assuming clash_value is a double
        double x = (static_cast<double>(coor[0]) + 0.5) * paras["grid_space"] + pocket[0][0];
        double y = (static_cast<double>(coor[1]) + 0.5) * paras["grid_space"] + pocket[1][0];
        double z = (static_cast<double>(coor[2]) + 0.5) * paras["grid_space"] + pocket[2][0];
        double vdw = 0;
        double elec = 0;

        for (const auto& atom : atoms) {
            double dist_square = pow(x - atom.coord[0], 2) + 
                                 pow(y - atom.coord[1], 2) + 
                                 pow(z - atom.coord[2], 2);
            
            if (dist_square <= paras["nb_cutoff_2"]) {
                elec += atom.charge / dist_square;
                double temp = paras["rmin_6"] / pow(dist_square, 3);
                vdw += paras["eps"] * temp * (temp - 2);
            }
        }
        return {vdw, elec};
    }
    return {value.vdw_or_clash, 0.0};
}

GRID Grid(const std::string& fpro, const VolRegion& pocket, int processes) {
    int eb_clash = static_cast<int>(paras["dist_clash"] / paras["grid_space"]);
    int ub_hbond = static_cast<int>(paras["ub_hbond_dist"] / paras["grid_space"]);
    double dist_clash_2 = pow(paras["dist_clash"], 2);
    double grid_space_2 = pow(paras["grid_space"], 2);

    std::array<std::array<double, 2>, 3> nbsize = pocket;
    nbsize[0][0] -= paras["nb_cutoff"];
    nbsize[0][1] += paras["nb_cutoff"];
    nbsize[1][0] -= paras["nb_cutoff"];
    nbsize[1][1] += paras["nb_cutoff"];
    nbsize[2][0] -= paras["nb_cutoff"];
    nbsize[2][1] += paras["nb_cutoff"];

    std::vector<TerAtom> atoms;
    std::array<int, 3> dim;
    dim[0] = static_cast<int>((pocket[0][1] - pocket[0][0]) / paras["grid_space"]) + 1;
    dim[1] = static_cast<int>((pocket[1][1] - pocket[1][0]) / paras["grid_space"]) + 1;
    dim[2] = static_cast<int>((pocket[2][1] - pocket[2][0]) / paras["grid_space"]) + 1;

    std::vector<std::vector<std::vector<GridValue>>> grid(
        dim[0],
        std::vector<std::vector<GridValue>>(
            dim[1],
            std::vector<GridValue>(
                dim[2],
                GridValue(0.0, 0.0, 0.0, 0.0)
            )
        )
    );
    
    std::ifstream f(fpro);
    std::string line;
    std::string atp ;
    while (std::getline(f, line)) {
        if (line.substr(0, 4) == "ATOM") {
            if (line[13] != 'H' && (line[16] == ' ' || line[16] == 'A')) {
                std::array<double, 3> coor = {
                    std::stod(line.substr(30, 8)),
                    std::stod(line.substr(38, 8)),
                    std::stod(line.substr(46, 8))
                };
                if (coor[0] >= nbsize[0][0] && coor[0] <= nbsize[0][1] &&
                    coor[1] > nbsize[1][0] && coor[1] <= nbsize[1][1] &&
                    coor[2] >= nbsize[2][0] && coor[2] <= nbsize[2][1]) {
                    
                    std::string res = line.substr(17, 3);
                    if (res == "HIS") res = "HSD";
                    atp = res + "_" + line.substr(13, 3);
                    double q = (atomtypes.count(atp) > 0) ? paras["elec_scaling"] * atomtypes[atp] : 0;
                    atoms.push_back({coor, q});
                }
                if (coor[0] >= pocket[0][0] && coor[0] <= pocket[0][1] &&
                    coor[1] > pocket[1][0] && coor[1] <= pocket[1][1] &&
                    coor[2] >= pocket[2][0] && coor[2] <= pocket[2][1]) {
                    
                    std::array<int, 3> loc = {
                        static_cast<int>((coor[0] - pocket[0][0]) / paras["grid_space"]),
                        static_cast<int>((coor[1] - pocket[1][0]) / paras["grid_space"]),
                        static_cast<int>((coor[2] - pocket[2][0]) / paras["grid_space"])
                    };

                    int xmin = std::max(loc[0] - eb_clash, 0);
                    int xmax = std::min(loc[0] + eb_clash, dim[0] - 1);
                    int ymin = std::max(loc[1] - eb_clash, 0);
                    int ymax = std::min(loc[1] + eb_clash, dim[1] - 1);
                    int zmin = std::max(loc[2] - eb_clash, 0);
                    int zmax = std::min(loc[2] + eb_clash, dim[2] - 1);

                    for (int i = xmin; i <= xmax; ++i) {
                        for (int j = ymin; j <= ymax; ++j) {
                            for (int k = zmin; k <= zmax; ++k) {
                                double dist = ((i - loc[0]) * (i - loc[0]) + 
                                               (j - loc[1]) * (j - loc[1]) + 
                                               (k - loc[2]) * (k - loc[2])) * grid_space_2;
                                if (dist <= dist_clash_2) {
                                    grid[i][j][k].vdw_or_clash = paras["e_clash"]; //记录clash能量，如果是0，就在CalcEn中计算vdw和elec
                                }
                            }
                        }
                    }

                    // Handle hydrogen bonds
                    if (hbondtypes.count(atp) > 0) {
                        xmin = std::max(loc[0] - ub_hbond, 0);
                        xmax = std::min(loc[0] + ub_hbond, dim[0] - 1);
                        ymin = std::max(loc[1] - ub_hbond, 0);
                        ymax = std::min(loc[1] + ub_hbond, dim[1] - 1);
                        zmin = std::max(loc[2] - ub_hbond, 0);
                        zmax = std::min(loc[2] + ub_hbond, dim[2] - 1);

                        for (int i = xmin; i <= xmax; ++i) {
                            for (int j = ymin; j <= ymax; ++j) {
                                for (int k = zmin; k <= zmax; ++k) {
                                    if (grid[i][j][k].vdw_or_clash == 0) { // 如果grid[i][j][k]的vdw_or_clash为0，则计算hbond
                                        double dist = std::sqrt((i - loc[0]) * (i - loc[0]) + 
                                                                 (j - loc[1]) * (j - loc[1]) + 
                                                                 (k - loc[2]) * (k - loc[2])) * paras["grid_space"];
                                        auto it = hbondtypes.find(atp);
                                        int hd_type;
                                        if (it != hbondtypes.end()) {
                                            hd_type = it->second;
                                            if (dist < paras["ub_hbond_dist"]) {
                                                if (hd_type == 2){ //如果蛋白原子是donor，则格点应该储存氢键受体的能量
                                                    if (grid[i][j][k].hd_acceptor.has_value()) {
                                                        grid[i][j][k].hd_acceptor.value() += (paras["ub_hbond_dist"] - dist) * paras["e_hbond"];
                                                        if (grid[i][j][k].hd_acceptor.value() < paras["e_hbond"]) {
                                                            grid[i][j][k].hd_acceptor.value() = paras["e_hbond"];
                                                        }
                                                    }
                                                }
                                                else if (hd_type == 3){ //如果蛋白原子是acceptor，则格点应该储存氢键给体的能量
                                                    if (grid[i][j][k].hd_donor.has_value()) {
                                                        grid[i][j][k].hd_donor.value() += (paras["ub_hbond_dist"] - dist) * paras["e_hbond"];
                                                        if (grid[i][j][k].hd_donor.value() < paras["e_hbond"]) {
                                                            grid[i][j][k].hd_donor.value() = paras["e_hbond"];
                                                        }
                                                    }
                                                }
                                                else {
                                                    std::cerr << "Error: hbondtype [" << hd_type << "] is not supported for atom [" << atp << "]" << std::endl;
                                                }
                                            }
                                        }
                                        else {
                                            grid[i][j][k].hd_acceptor = 0;
                                            grid[i][j][k].hd_donor = 0;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::mutex grid_mutex;
    std::vector<std::thread> threads;
    threads.reserve(processes);

    auto worker = [&](int start, int end) {
        for (int i = start; i < end; ++i) {
            for (int j = 0; j < dim[1]; ++j) {
                for (int k = 0; k < dim[2]; ++k) {
                    std::array<int, 3> coor = {i, j, k};
                    GridValue value = grid[i][j][k];
                    auto result = CalcEn(atoms, coor, pocket, value);
                    std::lock_guard<std::mutex> lock(grid_mutex);
                    grid[i][j][k] = GridValue{result[0], result[1], value.hd_donor, value.hd_acceptor};
                }
            }
        }
    };

    int chunk_size = dim[0] / processes;
    for (int i = 0; i < processes; ++i) {
        int start = i * chunk_size;
        int end = (i == processes - 1) ? dim[0] : (i + 1) * chunk_size;
        threads.emplace_back(worker, start, end);
    }

    for (auto& thread : threads) {
        thread.join();
    }

    return std::make_pair(grid, pocket);
}
