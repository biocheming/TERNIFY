#include <fstream>
#include <sstream>
#include <iostream>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <RDGeneral/Invariant.h>
#include <Geometry/point.h>

#include "ternify.hpp"


void read_parameters(const std::string& filename, Parameters& params) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open parameter file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;
        
        if (std::getline(iss, key, ':') && std::getline(iss, value)) {
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            if (key == "Interface") {
                std::istringstream value_stream(value);
                std::string number;
                std::vector<double> numbers;
                
                // 解析所有数字
                while (std::getline(value_stream, number, ',')) {
                    number.erase(0, number.find_first_not_of(" \t"));
                    number.erase(number.find_last_not_of(" \t") + 1);
                    if (!number.empty()) {
                        numbers.push_back(std::stod(number));
                    }
                }

                // 检查是否有足够的数字
                if (numbers.size() >= 6) {
                    for (int i = 0; i < 3; ++i) {
                        params.interface[i][0] = numbers[i * 2];
                        params.interface[i][1] = numbers[i * 2 + 1];
                    }
                } else {
                    throw std::runtime_error("Interface parameter requires 6 numbers");
                }
            }
            else if (key == "PROTACs") params.protacs_file = value;
            else if (key == "Warhead_anchor") params.warhead_anchor_file = value;
            else if (key == "Warhead_flex") params.warhead_flex_file = value;
            else if (key == "Protein_anchor") params.protein_anchor_file = value;
            else if (key == "Protein_flex") params.protein_flex_file = value;
            else if (key == "Output_protac") params.output_protac_file = value;
            else if (key == "Output_protein") params.output_protein_file = value;
            else if (key == "N_ini") params.n_ini = std::stoi(value);
            else if (key == "N_search") params.n_search = std::stoi(value);
            else if (key == "N_keep") params.n_keep = std::stoi(value);
            else if (key == "N_processes") params.n_processes = std::stoi(value);
        }
    }
}

void run_ternify(const Parameters& params) {
    // 读取warhead分子
    RDKit::SDMolSupplier w_anch_supplier(params.warhead_anchor_file,true,true,true); //读取的时候，分子并没有较好地检查，应该是bug
    RDKit::SDMolSupplier w_flex_supplier(params.warhead_flex_file,true,true, true);
    
    if (w_anch_supplier.atEnd() || w_flex_supplier.atEnd()) {
        throw std::runtime_error("Failed to read warhead molecules");
    }
    
    std::unique_ptr<RDKit::ROMol> w_anch(w_anch_supplier.next());
    std::unique_ptr<RDKit::ROMol> w_flex(w_flex_supplier.next());

    if (w_anch && w_flex) {
            std::cout << "READ: Number of atoms in WANCHOR: " << w_anch->getNumAtoms() << std::endl;
            std::cout << "READ: Number of bonds in WANCHOR: " << w_anch->getNumBonds() << std::endl;
            std::cout << "READ: Number of atoms in WFLEX: " << w_flex->getNumAtoms() << std::endl;
            std::cout << "READ: Number of bonds in WFLEX: " << w_flex->getNumBonds() << std::endl;
    }
    else {
        throw std::runtime_error("Failed to read warhead molecules");
    }
    // 创建输出文件
    RDKit::SDWriter protac_writer(params.output_protac_file);
    std::unique_ptr<std::ofstream> protein_writer;
    bool write_protein = params.output_protein_file != "None";
    
    if (write_protein) {
        protein_writer.reset(new std::ofstream(params.output_protein_file));
    }

    // 创建网格
    std::cout << "Creating grid for anchor_pro..." << std::endl;
    GRID grid_anchor;
    bool grid_anchor_initialized = false;
    try {
        grid_anchor = Grid(params.protein_anchor_file, params.interface, params.n_processes);
        grid_anchor_initialized = true;
    } catch (const std::exception& e) {
        std::cerr << "Failed to create grid for anchor_pro: " << e.what() << std::endl;
    }
    // 检查网格是否已初始化
    if (!grid_anchor_initialized) {
        std::cerr << "Warning: grid_anchor was not initialized." << std::endl;
    }
    
    // 计算柔性蛋白网格的位置
    std::cout << "Creating grid for flex_pro..." << std::endl;
    const RDKit::Conformer& conf_flex = w_flex->getConformer();
    VolRegion site_flex;
    
    for (int i = 0; i < 3; ++i) {
        double min_val = std::numeric_limits<double>::max();
        double max_val = std::numeric_limits<double>::lowest();
        
        for (unsigned int j = 0; j < conf_flex.getNumAtoms(); ++j) {
            const RDGeom::Point3D& pos = conf_flex.getAtomPos(j);
            double coord = (i == 0) ? pos.x : ((i == 1) ? pos.y : pos.z);
            min_val = std::min(min_val, coord);
            max_val = std::max(max_val, coord);
        }
        
        site_flex[i][0] = min_val - 5.0;  // 使用数组索引而不是 push_back
        site_flex[i][1] = max_val + 5.0;    
    }

    GRID grid_flex;
    bool grid_flex_initialized = false;
    try {
        grid_flex = Grid(params.protein_flex_file, site_flex, params.n_processes);
        grid_flex_initialized = true;
    } catch (const std::exception& e) {
        std::cerr << "Failed to create grid for flex_pro: " << e.what() << std::endl;
    }


    if (!grid_flex_initialized) {
        std::cerr << "Warning: grid_flex was not initialized." << std::endl;
    }
    // 处理PROTAC分子
    std::cout << "Reading protac file..." << std::endl;
    RDKit::SDMolSupplier protacs_supplier(params.protacs_file,true,true,false);
    size_t total_protacs = protacs_supplier.length();
    std::cout << "There is(are) " << protacs_supplier.length() << " protac molecule(s) to be processed..." << std::endl;
    std::unique_ptr<RDKit::ROMol> mol;
    for (size_t protac_index = 0; protac_index < total_protacs; ++protac_index) {
        mol = std::unique_ptr<RDKit::ROMol>(protacs_supplier[protac_index]);
        // Check if the molecule was successfully read and matches the substructure
        if (mol) {
            std::cout << "READ: Number of atoms in protac: " << mol->getNumAtoms() << std::endl;
            std::cout << "READ: Number of bonds in protac: " << mol->getNumBonds() << std::endl;
            if (RDKit::SubstructMatch(*mol, *w_anch).empty() || 
                RDKit::SubstructMatch(*mol, *w_flex).empty() ){
                std::cout << "WARNING: Skipping molecule " << protac_index + 1 << ": this protac does not match substructures with anch/flex." << std::endl;
                continue;
            }
            Protac PROTac;
            std::cout << "Initializing protac..." << std::endl;
            // Initialize Protac object with the current molecule and other parameters
            PROTac.init(mol.get(), w_anch.get(), w_flex.get(),
                    params.protein_flex_file, params.n_processes,
                    grid_anchor, grid_flex);

            // Print Protac information
            //PROTac.printProtacInfo();

            std::cout << "Sampling and Searching..." << std::endl;
            PROTac.sample(params.n_ini, params.n_search);
        
            std::cout << "Writing output..." << std::endl;
            // Write the output to the specified writer
            PROTac.output(protac_writer, protein_writer ? *protein_writer : std::cout, params.n_keep, write_protein);            
        } 
        else {
            std::cout << "Warning: Molecule " << protac_index + 1 << " is invalid or could not be read." << std::endl;
        }
        
    }

}

int main(int argc, char* argv[]) {
    try {
        // 命令行参数解析
        std::string param_file = "tcs.inp";
        if (argc > 1) {
            for (int i = 1; i < argc; i++) {
                std::string arg = argv[i];
                if (arg == "-p" && i + 1 < argc) {
                    param_file = argv[++i];
                }
            }
        }

        // 打印版本信息
        std::cout << "TERNIFY: Efficient Sampling of PROTAC-Induced Ternary Complexes\n"
                  << "Hongtao Zhao, PhD\n"
                  << "Ximing XU, PhD [C++ implementation]\n"
                  << "Version: 2024-11-17" << std::endl;

        // 读取参数并运行
        Parameters params;
        std::cout << "Reading parameters from " << param_file << std::endl;
        read_parameters(param_file, params);
        run_ternify(params);

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}