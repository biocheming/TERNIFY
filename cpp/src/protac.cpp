#include <vector>
#include <memory>
#include <array>
#include <iterator>
#include <iomanip>  

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>

#include <RDGeneral/types.h>     // for RDGeom::POINT3D_VECT
#include <Geometry/point.h>      // for RDGeom::Point3D

#include <Eigen/Dense>
#include "protac.hpp"
#include "minimize_h.hpp"
#include "lowEnerDihe.hpp"
// 定义thread_local变量
thread_local std::unique_ptr<RDKit::ROMol> Protac::thread_local_mol_;
std::mutex Protac::output_mutex_;

void Protac::init(RDKit::ROMol* protac,  // unaligned protac
                 RDKit::ROMol* w_anch,   //E3 lig
                 RDKit::ROMol* w_flex,   //POI lig
                 const std::string& fpro_flex, // POI protein
                 int processes,
                 const GRID& grid_anchor,
                 const GRID& grid_flex,
                 bool verbose) {
    
    processes_ = processes;
    grid_anchor_ = grid_anchor;
    grid_flex_ = grid_flex;

    // Initialize random number generators
    rng_ = std::mt19937(std::random_device{}());
    angle_dist_ = std::uniform_real_distribution<double>(-180.0, 180.0);
    uniform_dist_ = std::uniform_real_distribution<double>(0.0, 1.0);

    // 初始化分子, 并优化
    protac_= MinimizeH(*protac, 1.0e4, false);
    RDKit::RWMol protac_ini_H = *protac_ ;
    
    // 获取构象
    RDKit::Conformer& conf_protac = protac_->getConformer();        //protac whole part
    const RDKit::Conformer& conf_anch = w_anch->getConformer(); 
    const RDKit::Conformer& conf_flex = w_flex->getConformer();

    // 移动w_flex坐标到原点
    RDGeom::POINT3D_VECT coords_flex_conf = conf_flex.getPositions();
    RDGeom::Point3D translation(0.0, 0.0, 0.0);
    for (const auto& point : coords_flex_conf) {
        translation += point;
    }
    translation /= coords_flex_conf.size();
    translation_ = {translation.x, translation.y, translation.z};  // 保存平移向量
    // 应用平移
    for (auto& point : coords_flex_conf) {
        point -= translation;
        coord_subs_var.push_back({point.x, point.y, point.z});
    }
    // 读取蛋白质，传入平移向量, w_flex, fpro_flex都不带H
    protein_.ReadProt(fpro_flex, translation_); 
    // 对齐柔性warhead
    // 设置更宽松的匹配参数
    RDKit::SubstructMatchParameters match_params;
    match_params.useChirality = false;       // 不考虑手性
    match_params.useQueryQueryMatches = true; // 允许查询-查询匹配
    match_params.maxMatches = 1000;          // 允许多个匹配
    match_params.uniquify = true;            // 确保匹配是唯一的
    match_params.recursionPossible = true;   // 允许递归匹配
    match_params.aromaticMatchesConjugated = true; // 芳香族可以匹配共轭
    

    std::vector<RDKit::MatchVectType> matches_01;
    matches_01 = RDKit::SubstructMatch(*protac_, *w_flex, match_params);        
    if (matches_01.empty()) {
        throw std::runtime_error("No match found for the flexible warhead.");
    }
    RDKit::MatchVectType match_01 = matches_01[0];
    
    // 保存protac_弹头原子索引用于RMSD计算
    warhead_atoms_.clear();
    for (const auto& pair : match_01) {
        warhead_atoms_.push_back(pair.second);  // 保存protac中弹头原子的索引
    }
    
    // 查找连接点
    int cp_01 = -1;
    int cp_01_neighbor = -1;  // 添加这行来存储实际的连接原子
    for (size_t i = 0; i < match_01.size(); ++i) {
        RDKit::Atom* atom = protac_->getAtomWithIdx(match_01[i].second);  
        // 获取邻接原子的迭代器对
        auto [nbrIdx, endNbrs] = protac_->getAtomNeighbors(atom);
        // 遍历所有邻居
        while (nbrIdx != endNbrs) {
            int neighborIdx = *nbrIdx;
            // 检查邻居是否不在 match_01 中
            if (std::find_if(match_01.begin(), match_01.end(),
                [&](const auto& p) { return p.second == neighborIdx; }) 
                    == match_01.end() && protac_->getAtomWithIdx(neighborIdx)->getAtomicNum() != 1) {
                cp_01 = match_01[i].second;  // 存储实际的原子索引
                cp_01_neighbor = neighborIdx;  // 存储邻居原子的索引
                break;
            }
            ++nbrIdx;
        }
        if (cp_01 != -1) break;
    }
    
    if (cp_01 == -1) {
        std::cerr << "Could not find a connected unmatched atom." << std::endl;
        return;
    } else {
        if (verbose) {
            std::cout << "Connection 1: [" << cp_01 << "] - [" << cp_01_neighbor << "]" << std::endl;
        }
    }

    // 检查匹配的有效性
    for (const auto& pair : match_01) {
        if (static_cast<unsigned int>(pair.first) >= w_flex->getNumAtoms() ||
            static_cast<unsigned int>(pair.second) >= protac_->getNumAtoms()) {
            throw std::runtime_error("Invalid atom indices in match: " + 
                                   std::to_string(pair.first) + ", " + 
                                   std::to_string(pair.second));
        }
    }
    // 新的对齐策略：使用最小匹配单元进行精准对齐
    if (verbose) {
        std::cout << "Creating minimal matching unit for precise alignment..." << std::endl;
    }
    
    // 找到连接点cp_01在w_flex中对应的原子索引
    int cp_01_in_flex = -1;
    for (const auto& pair : match_01) {
        if (pair.second == cp_01) {
            cp_01_in_flex = pair.first;
            break;
        }
    }
    
    if (cp_01_in_flex == -1) {
        throw std::runtime_error("Could not find cp_01 in w_flex match");
    }
    
    // 构建最小匹配单元：cp_01 + 与其在弹头中有键连接的原子
    std::vector<std::pair<int, int>> minimal_pairs;
    minimal_pairs.push_back(std::make_pair(cp_01, cp_01_in_flex));  // 连接点本身
    
    // 找到与cp_01在w_flex中有键连接的原子
    for (const auto& pair : match_01) {
        int protac_idx = pair.second;
        int flex_idx = pair.first;
        
        // 检查w_flex中是否存在cp_01_in_flex与flex_idx之间的键
        if (flex_idx != cp_01_in_flex && 
            w_flex->getBondBetweenAtoms(cp_01_in_flex, flex_idx)) {
            minimal_pairs.push_back(std::make_pair(protac_idx, flex_idx));
        }
    }

    if (verbose) {
        std::cout << "Minimal matching unit size: " << minimal_pairs.size() << " atoms" << std::endl;
        for (const auto& pair : minimal_pairs) {
            std::cout << "  protac atom " << pair.first << " -> w_flex atom " << pair.second << std::endl;
        }
    }
    
    // 创建对齐匹配向量（注意顺序：protac -> w_flex）
    RDKit::MatchVectType aligned_match;
    for (const auto& pair : minimal_pairs) {
        aligned_match.push_back(std::make_pair(pair.first, pair.second));
    }

    // protac_是probe（要移动的），w_flex是ref（固定的参考）
    std::cout << "Aligning protac to flexible warhead..." << std::endl;
    try {
        // 使用最小匹配单元进行对齐
        RDKit::MolAlign::alignMol(*protac_, *w_flex, -1, -1, &aligned_match);
        std::cout << "Alignment completed" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error during alignment: " << e.what() << std::endl;
        throw;
    }
    // 设置原子位置并记录索引
    // 首先记录要移动的重原子的原始坐标
    std::vector<RDGeom::Point3D> original_positions;
    for (size_t i = 0; i < match_01.size(); ++i) {
        original_positions.push_back(conf_protac.getAtomPos(match_01[i].second));
    }
    // 更新原子位置
    for (size_t i = 0; i < match_01.size(); ++i) {
        const RDGeom::Point3D& coor = conf_flex.getAtomPos(match_01[i].first); // 使用match_01[i].first来获取w_flex中的原子坐标
        const RDGeom::Point3D& original_heavy_pos = original_positions[i]; 

        conf_protac.setAtomPos(match_01[i].second, coor);
        idx_.push_back(match_01[i].second);
        
        // 计算平移向量并应用到连接的氢原子
        RDGeom::Point3D translation = coor - original_heavy_pos;
        RDKit::Atom* heavy_atom = protac_->getAtomWithIdx(match_01[i].second);
        
        // 获取邻接原子的迭代器
        auto [nbrIdx, endNbrs] = protac_->getAtomNeighbors(heavy_atom);
        while (nbrIdx != endNbrs) {
            RDKit::Atom* neighbor = protac_->getAtomWithIdx(*nbrIdx);
            if (neighbor->getAtomicNum() == 1) {  // 氢原子
                RDGeom::Point3D h_pos_before = conf_protac.getAtomPos(*nbrIdx);
                RDGeom::Point3D h_target_pos = h_pos_before + translation;
                conf_protac.setAtomPos(*nbrIdx, h_target_pos);
            }
            ++nbrIdx;
        }
    }
    //optimizeWithFixedAtoms(*protac_, idx_); //protac_的第2次优化
    //optimizeWithConstrAtoms(*protac_, idx_);
    optimizeH(*protac_, 1.0e3);
    conf_protac = protac_->getConformer();
    std::cout << "Finding acceptors..." << std::endl;
    // 识别氢键受体
    std::unique_ptr<RDKit::ROMol> acceptor_patt(RDKit::SmartsToMol("[#8,nX2,$(N#C);+0]"));
    std::vector<int> hb_acceptors;
    std::vector<RDKit::MatchVectType> hb_acc_matches;
    RDKit::SubstructMatch(*protac_, *acceptor_patt, hb_acc_matches);
    for (const auto& acc_match : hb_acc_matches) {
        for (const auto& pair : acc_match) {
            hb_acceptors.push_back(pair.second);
        }
    }

    std::cout << "Finding donors..." << std::endl;
    // 识别氢键供体
    std::unique_ptr<RDKit::ROMol> donor_patt(RDKit::SmartsToMol("[#7,#8;H;+0]"));
    std::vector<int> hb_donors;
    std::vector<RDKit::MatchVectType> hb_don_matches;
    RDKit::SubstructMatch(*protac_, *donor_patt, hb_don_matches);
    for (const auto& don_match : hb_don_matches) {
        for (const auto& pair : don_match) {
            hb_donors.push_back(pair.second);
        }
    }

    // 先进行第一次对齐，但暂不计算电荷
    // 计算不在match_01中的原子的索引，稍后计算电荷
    std::vector<int> match_01_indices;
    for (const auto& pair : match_01) {
        match_01_indices.push_back(pair.second);
    }
    
    // 对齐锚定warhead w_anch
    std::vector<RDKit::MatchVectType> matches_02;
    matches_02 = RDKit::SubstructMatch(*protac_, *w_anch, match_params);
    if(matches_02.empty()){
        throw std::runtime_error("No match found for the anchor warhead.");
    }
    RDKit::MatchVectType match_02 = matches_02[0];

    int cp_02 = -1;
    int cp_02_neighbor = -1;
    for (size_t i = 0; i < match_02.size(); ++i) {
        RDKit::Atom* atom = protac_->getAtomWithIdx(match_02[i].second);
        // 获取邻接原子的迭代器对
        auto [nbrIdx, endNbrs] = protac_->getAtomNeighbors(atom);
        // 遍历所有邻居
        while (nbrIdx != endNbrs) {
            int neighborIdx = *nbrIdx;
            // 检查邻居是否不在 match_02 中
            if (std::find_if(match_02.begin(), match_02.end(),
                [&](const auto& p) { return p.second == neighborIdx; }) 
                == match_02.end() && protac_->getAtomWithIdx(neighborIdx)->getAtomicNum() != 1) {
                cp_02 = match_02[i].second;  // 存储实际的原子索引
                cp_02_neighbor = neighborIdx;  // 存储邻居原子的索引
                break;
            }
            ++nbrIdx;
        }
        if (cp_02 != -1) break;
    }

    // 添加调试信息
    if (cp_02 != -1) {
        if (verbose) {
            std::cout << "Connection 2: [" << cp_02 << "] - [" << cp_02_neighbor << "]" << std::endl;
        }
    } else {
        std::cerr << "Could not find a connected unmatched atom." << std::endl;
        return;
    }
    

   // 新的对齐策略：使用最小匹配单元进行精准对齐（第二个弹头）
    std::cout << "Creating minimal matching unit for anchor alignment..." << std::endl;
    
    if (cp_02 == -1) {
        throw std::runtime_error("Could not find connection point cp_02");
    }
    
    // 构建最小匹配单元：按照Python版本的逻辑
    std::vector<std::pair<int, int>> minimal_pairs_anch;
    
    // 首先添加连接点：(protac中的原子索引, w_anch中的原子索引)
    // 这里需要找到cp_02在match_02中的位置
    int cp_02_match_index = -1;
    for (size_t i = 0; i < match_02.size(); ++i) {
        if (match_02[i].second == cp_02) {  // cp_02是protac中的原子索引
            cp_02_match_index = i;
            break;
        }
    }
    
    if (cp_02_match_index == -1) {
        throw std::runtime_error("Could not find cp_02 in match_02");
    }
    
    // 添加连接点
    minimal_pairs_anch.push_back(std::make_pair(
        match_02[cp_02_match_index].second,  // protac中的原子索引
        match_02[cp_02_match_index].first    // w_anch中的原子索引
    ));
    if (verbose) {
        std::cout << "Connection point: protac atom " << match_02[cp_02_match_index].second 
                << " -> w_anch atom " << match_02[cp_02_match_index].first << std::endl;
    }
    
    // 找到与连接点在w_anch中有键连接的原子
    for (size_t i = 0; i < match_02.size(); ++i) {
        if (static_cast<int>(i) != cp_02_match_index) {  // 不是连接点本身
            // 检查w_anch中是否存在连接点与当前原子之间的键
            if (w_anch->getBondBetweenAtoms(match_02[cp_02_match_index].first, match_02[i].first)) {
                minimal_pairs_anch.push_back(std::make_pair(
                    match_02[i].second,  // protac中的原子索引
                    match_02[i].first    // w_anch中的原子索引
                ));
                if (verbose) {  
                    std::cout << "Found bonded atom: protac atom " << match_02[i].second 
                            << " -> w_anch atom " << match_02[i].first << std::endl;
                }
            }
        }
    }
    
    // 创建对齐匹配向量
    aligned_match.clear();
    for (const auto& pair : minimal_pairs_anch) {
        aligned_match.push_back(std::make_pair(pair.first, pair.second));
    }
    
    // 执行分子对齐
    std::cout << "Aligning protac to anchor warhead using minimal matching unit..." << std::endl;
    try {
        RDKit::MolAlign::alignMol(*protac_, *w_anch, -1, -1, &aligned_match);
        std::cout << "Anchor minimal unit alignment completed" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during anchor minimal unit alignment: " << e.what() << std::endl;
        throw;
    }

    // 更新原子位置
    std::vector<int> AnchMatchedIDInProtac;
    // 首先记录要移动的重原子的原始坐标
    std::vector<RDGeom::Point3D> original_anch_positions;
    for (size_t i = 0; i < match_02.size(); ++i) {
        original_anch_positions.push_back(conf_protac.getAtomPos(match_02[i].second));
    }
        
    for (size_t i = 0; i < match_02.size(); ++i) {
        const RDGeom::Point3D& coor = conf_anch.getAtomPos(match_02[i].first);
        // Original heavy atom position in protac_ before moving IT
        const RDGeom::Point3D& original_heavy_pos = original_anch_positions[i]; 
        
        conf_protac.setAtomPos(match_02[i].second, coor);
        AnchMatchedIDInProtac.push_back(match_02[i].second);

        RDGeom::Point3D translation = coor - original_heavy_pos;
        RDKit::Atom* heavy_atom = protac_->getAtomWithIdx(match_02[i].second);
        // 获取邻接原子的迭代器
        auto [nbrIdx, endNbrs] = protac_->getAtomNeighbors(heavy_atom);
        while (nbrIdx != endNbrs) {
            RDKit::Atom* neighbor = protac_->getAtomWithIdx(*nbrIdx);
            if (neighbor->getAtomicNum() == 1) {  // 氢原子
                RDGeom::Point3D h_pos_before = conf_protac.getAtomPos(*nbrIdx);
                RDGeom::Point3D h_target_pos = h_pos_before + translation;
                conf_protac.setAtomPos(*nbrIdx, h_target_pos);
            }
            ++nbrIdx;
        }
    }
    //MiniFixAtomTor(*protac_, AnchMatchedIDInProtac, idx_); 
    optimizeH(*protac_, 100.0);
    conf_protac = protac_->getConformer();
    
    // 所有分子对齐完成后重新计算Gasteiger电荷（匹配Python版本的时机）
//    RDKit::computeGasteigerCharges(*protac_);
    
    // 首先计算q_anchor：不在match_01中的原子（只包含重原子）
    for (size_t i = 0; i < protac_->getNumAtoms(); ++i) {
        if (std::find(match_01_indices.begin(), match_01_indices.end(), i) 
            == match_01_indices.end()) {
            RDKit::Atom* atom = protac_->getAtomWithIdx(i);
            // 跳过氢原子，只处理重原子（匹配Python版本）
            if (atom->getAtomicNum() == 1) continue;
            
            // 获取Gasteiger电荷
            double gasteiger_charge = 0.0;
            if (atom->hasProp("_GasteigerCharge")) {
                atom->getProp("_GasteigerCharge", gasteiger_charge);
            }
            
            // 手动累加连接氢原子的电荷
            double hydrogen_charges = 0.0;
            auto [nbrIdx, endNbrs] = protac_->getAtomNeighbors(atom);
            while (nbrIdx != endNbrs) {
                RDKit::Atom* neighbor = protac_->getAtomWithIdx(*nbrIdx);
                if (neighbor->getAtomicNum() == 1) {  // 氢原子
                    double h_charge = 0.0;
                    if (neighbor->hasProp("_GasteigerCharge")) {
                        neighbor->getProp("_GasteigerCharge", h_charge);
                        hydrogen_charges += h_charge;
                    }
                }
                ++nbrIdx;
            }
            
            double total_charge = gasteiger_charge + hydrogen_charges;
            // 确定氢键类型 - 使用前面计算的氢键信息
            std::optional<int> hb_type = std::nullopt;
            if (std::find(hb_donors.begin(), hb_donors.end(), i) != hb_donors.end()) {
                hb_type = 3;  // 氢键供体设为3（与protein.cpp一致）
            } else if (std::find(hb_acceptors.begin(), hb_acceptors.end(), i) != hb_acceptors.end()) {
                hb_type = 2;  // 氢键受体设为2（与protein.cpp一致）
            }
            
            // 添加到q_anchor列表
            hb_type = std::nullopt;
            q_anchor_.push_back(std::make_tuple(i, total_charge, hb_type));
        }
    }
    
    // 然后计算q_flex：不在match_02中的原子（只包含重原子）
    std::vector<int> match_02_indices;
    for (const auto& pair : match_02) {
        match_02_indices.push_back(pair.second);
    }
    
    for (size_t i = 0; i < protac_->getNumAtoms(); ++i) {
        if (std::find(match_02_indices.begin(), match_02_indices.end(), i) 
            == match_02_indices.end()) {
            RDKit::Atom* atom = protac_->getAtomWithIdx(i);
            // 跳过氢原子，只处理重原子（匹配Python版本）
            if (atom->getAtomicNum() == 1) continue;
            
            // 获取Gasteiger电荷
            double gasteiger_charge = 0.0;
            if (atom->hasProp("_GasteigerCharge")) {
                atom->getProp("_GasteigerCharge", gasteiger_charge);
            }
            
            // 手动累加连接氢原子的电荷
            double hydrogen_charges = 0.0;
            auto [nbrIdx, endNbrs] = protac_->getAtomNeighbors(atom);
            while (nbrIdx != endNbrs) {
                RDKit::Atom* neighbor = protac_->getAtomWithIdx(*nbrIdx);
                if (neighbor->getAtomicNum() == 1) {  // 氢原子
                    double h_charge = 0.0;
                    if (neighbor->hasProp("_GasteigerCharge")) {
                        neighbor->getProp("_GasteigerCharge", h_charge);
                        hydrogen_charges += h_charge;
                    }
                }
                ++nbrIdx;
            }
            
            double total_charge = gasteiger_charge + hydrogen_charges;
            // 确定氢键类型 - 使用前面计算的氢键信息
            std::optional<int> hb_type = std::nullopt;
            if (std::find(hb_donors.begin(), hb_donors.end(), i) != hb_donors.end()) {
                hb_type = 3;  // 氢键供体
            } else if (std::find(hb_acceptors.begin(), hb_acceptors.end(), i) != hb_acceptors.end()) {
                hb_type = 2;  // 氢键受体
            }
            
            // 添加到q_flex列表
            hb_type = std::nullopt;
            q_flex_.push_back(std::make_tuple(i, total_charge, hb_type));
        }
    }

    // 查找linker
    std::vector<int> linker;
    std::vector<int> warheads;
    
    // 合并match_01和match_02的索引
    for (const auto& pair : match_01) {
        warheads.push_back(pair.second);
    }
    for (const auto& pair : match_02) {
        warheads.push_back(pair.second);
    }
    
    // 扩展warheads，包括直连的氢原子
    std::vector<int> complete_warheads = warheads;  // 先复制重原子列表
    for (int warhead_atom : warheads) {
        RDKit::Atom* atom = protac_->getAtomWithIdx(warhead_atom);
        auto [nbrIdx, endNbrs] = protac_->getAtomNeighbors(atom);
        while (nbrIdx != endNbrs) {
            RDKit::Atom* neighbor = protac_->getAtomWithIdx(*nbrIdx);
            if (neighbor->getAtomicNum() == 1) {  // 氢原子
                complete_warheads.push_back(neighbor->getIdx());
            }
            ++nbrIdx;
        }
    }
    warheads = complete_warheads;

    // 找出不在warheads中的原子
    for (size_t i = 0; i < protac_->getNumAtoms(); ++i) {
        if (std::find(warheads.begin(), warheads.end(), i) == warheads.end()) {
            linker.push_back(i);
        }
    }

    // 查找可旋转键创建SMARTS模式
    std::cout << "Finding rotatable bonds..." << std::endl;
    std::string RotSmarts = "[!$(C(=[N,O,S])-!@[#7H,O,S])&!$([#7H,O,S]-!@C=[N,O,S])&!D1]-&!@[!D1]";
    RDKit::ROMol* RotPatt = RDKit::SmartsToMol(RotSmarts);

    if (!RotPatt) {
        throw std::runtime_error("Failed to parse SMARTS pattern");
    }

    // 查找匹配的键
    std::vector<RDKit::MatchVectType> RotMatches = RDKit::SubstructMatch(*protac_, *RotPatt);
    // 清理SMARTS模式
    delete RotPatt;
    RDKit::MatchVectType rbond;
    // 处理找到的键
    for (const auto& rot_match : RotMatches) {
        if (rot_match.size() >= 2) {
            int a = rot_match[0].second;
            int b = rot_match[1].second;
            
            if (std::find(linker.begin(), linker.end(), a) != linker.end() || 
                std::find(linker.begin(), linker.end(), b) != linker.end()) {
                
                // 直接使用cp_02，因为它已经是实际的原子索引
                int d1 = (a != cp_02) ? 
                         RDKit::MolOps::getShortestPath(*protac_, a, cp_02).size() : 0;

                int d2 = (b != cp_02) ? 
                         RDKit::MolOps::getShortestPath(*protac_, b, cp_02).size() : 0;

                if (d1 < d2) {
                    rbond.push_back(std::make_pair(a, b));
                } else {
                    rbond.push_back(std::make_pair(b, a));
                }
            }
        }
    }
    std::cout << "Finding dihedrals..." << std::endl;
    // 查找二面角
    for (const auto& bond : rbond) {
        std::vector<int> na;
        // 查找第一个原子的邻居
        RDKit::Atom* atom1 = protac_->getAtomWithIdx(bond.first);
        auto [nbrIdx, endNbrs] = protac_->getAtomNeighbors(atom1);
        // 遍历所有邻居
        while (nbrIdx != endNbrs) {
            int neighborIdx = *nbrIdx;
            if (neighborIdx != bond.second) {  // 修改：检查不是键的另一端
               na.push_back(neighborIdx);
               break;
            }
            ++nbrIdx;
        }
        // 查找第二个原子的邻居
        RDKit::Atom* atom2 = protac_->getAtomWithIdx(bond.second);
        auto [neighborIdx, endNeighborIdx] = protac_->getAtomNeighbors(atom2);

        // 遍历所有邻居
        while (neighborIdx != endNeighborIdx) {
            int neighborIndex = *neighborIdx;
            if (neighborIndex != bond.first) {  // 修改：检查不是键的另一端
                na.push_back(neighborIndex);
                break;
            }
            ++neighborIdx;
        }
    
        // 确保找到了两个邻居原子并验证它们都是不同的
        if (na.size() == 2 && na[0] != na[1] && 
            na[0] != bond.first && na[0] != bond.second && 
            na[1] != bond.first && na[1] != bond.second) {
            rot_dihe_.push_back({na[0], bond.first, bond.second, na[1]});
        
            } else {
                std::cerr << "Warning: Skipping invalid dihedral for bond " 
                            << bond.first << "-" << bond.second << std::endl;
                continue;
            }
    }                   
    // 调用list函数处理warheads和rbond
    std::cout << "Processing FF parameters..." << std::endl;
    list(warheads, rbond, verbose);
    //protac_ = MinimizeH(*protac_, 15.0, false);
    // 计算参考内部能量（基于经过两次对齐和优化的构象）
    E_intra_ref_ = e_intra(protac_.get());
    std::cout << "E_intra_ref: " << E_intra_ref_ << std::endl;
    // I want to check the linker atoms 
    if(verbose > 0 ){
        std::cout << "Linker atoms ( " << linker.size() << " ): ";
        for (int atom_idx : linker) {
            RDKit::Atom* atom_i = protac_->getAtomWithIdx(atom_idx);
            std::cout << atom_idx - 1 << "-" << atom_i->getSymbol() << " ";
        }
        std::cout << std::endl;
    }  
}

// ==================扫描单个二面角的低能构象=========================
std::vector<double> Protac::scanTorsion(double v1, double v2, double v3) {
    std::vector<double> lowEnerDiheAngles;
    std::vector<std::pair<double, double>> angle_energies;
    // 从-180°到180°，每隔5°扫描
    for (double angle = -180.0; angle <= 180.0; angle += 5.0) {
        double angle_rad = angle * M_PI / 180.0;
        
        // 计算该角度下的二面角能量
        double energy = 0.5 * (
            v1 * (1.0 + std::cos(angle_rad)) +
            v2 * (1.0 - std::cos(2.0 * angle_rad)) +
            v3 * (1.0 + std::cos(3.0 * angle_rad))
        );
        
        angle_energies.emplace_back(angle, energy);
    }
    
    // 按能量排序
    std::sort(angle_energies.begin(), angle_energies.end(),
             [](const auto& a, const auto& b) { return a.second < b.second; });
    
    // 只保留前36个最低能量对应的角度
    for (size_t i = 0; i < std::min(size_t(36), angle_energies.size()); ++i) {
        lowEnerDiheAngles.push_back(angle_energies[i].first);
    }
        
    return lowEnerDiheAngles;
}

std::vector<double> Protac::listTorsion(double grid_step){
    std::vector<double> listDiheAngles;
    for (double angle = -180.0; angle <= 180.0; angle += grid_step) {
        listDiheAngles.push_back(angle);
    }    
    return listDiheAngles;
}

//=================== GET FF PARAMETERS ========================
void Protac::list(const std::vector<int>& warheads, const std::vector<std::pair<int, int>>& rbond, bool print_info) {
    // 首先检查输入参数的有效性
    if (!protac_) {
        throw std::runtime_error("Protac molecule is null");
    }

    // 检查warheads和rbond是否为空
    if (warheads.empty()) {
        throw std::runtime_error("Warning: warheads vector is empty");
    }    
    // 创建warhead和氢原子的映射
    std::unordered_map<int, bool> wh_h;
    // 添加warhead原子和它们的氢原子
    for (unsigned int i : warheads) {
        if (i >= protac_->getNumAtoms()) {
            std::cerr << "Invalid atom index: " << i << std::endl;
            continue;
        }

        wh_h[i] = true;
        RDKit::Atom* atom = protac_->getAtomWithIdx(i);

        if (atom) {
            // 获取邻接原子的迭代器
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = protac_->getAtomNeighbors(atom);

            // 遍历邻居原子
            while (nbrIdx != endNbrs) {
                RDKit::Atom* neighbor = protac_->getAtomWithIdx(*nbrIdx);
                if (neighbor->getAtomicNum() == 1) {  // 检查邻居是否是氢原子
                    wh_h[neighbor->getIdx()] = true;
                }
                ++nbrIdx;  // 移动到下一个邻居
            }
        }
    }

    // 获取MMFF属性
    RDKit::MMFF::MMFFMolProperties mmffProps(*protac_, "MMFF94s");
    if (!mmffProps.isValid()) {
        throw std::runtime_error("Failed to setup MMFF properties");
    }

    // 获取拓扑距离矩阵，而不是3D距离矩阵
    const double* topDistMat = RDKit::MolOps::getDistanceMat(*protac_);
    
    size_t n = protac_->getNumAtoms();
    std::vector<std::vector<double>> topDist(n, std::vector<double>(n));  // 创建一个二维向量来存储拓扑距离矩阵

    // 将指针中的数据复制到二维向量中
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            topDist[i][j] = topDistMat[i * n + j];  // 计算索引并赋值
        }
    }

    // 打印拓扑距离矩阵（可选）
    if(print_info){
        std::cout << "Topological distance matrix (first 5x5):" << std::endl;
        for (size_t i = 0; i < std::min(n, size_t(5)); ++i) {
            for (size_t j = 0; j < std::min(n, size_t(5)); ++j) {
                std::cout << topDist[i][j] << "\t";
            }
            std::cout << std::endl;
        }
    }

    // 计算 VdW 参数
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            // 只考虑拓扑距离 >= 3 的原子对（非共价相互作用）
            if (topDist[i][j] > 2.0) {  
                // 如果两个原子都在wh_h中，跳过
                if (wh_h.find(i) != wh_h.end() && wh_h.find(j) != wh_h.end()) {
                    continue;
                }

                // 获取VdW参数
                RDKit::MMFF::MMFFVdWRijstarEps vdwParams;
                if (!mmffProps.getMMFFVdWParams(i, j, vdwParams)) {
                    continue;
                }
                
                // 创建VdW参数列表
                VdwParam param;
                param.atom1_idx = i;
                param.atom2_idx = j;
                param.param1 = vdwParams.epsilon;
                param.param1_1 = vdwParams.R_ij_star;
                param.param2 = 1.07 * vdwParams.R_ij_star;
                param.param3 = 0.07 * vdwParams.R_ij_star;
                double r_ij_star_pow7 = std::pow(vdwParams.R_ij_star, 7);  // 先计算7次方
                param.param4 = 1.12 * r_ij_star_pow7;  // 再乘以系数
                param.param5 = 0.12 * r_ij_star_pow7;  // 再乘以系数
                
                // 对拓扑距离=3的1-4连接键进行能量缩放（一般是0.75）
                if (std::abs(topDist[i][j] - 3.0) < 0.001) {  // 使用浮点比较
                    param.param1 *= 0.75;  // 对epsilon进行缩放
                }
                
                // 打印VdW参数
                /*
                if (print_info) {
                    
                    std::cout << std::fixed << std::setprecision(6)
                              << "VdW_params: atoms=" << i << "-" << j 
                              << " epsilon=" << param.param1 
                              << " R_ij_star=" << param.param1_1
                              << " param2=" << param.param2
                              << " param3=" << param.param3
                              << " param4=" << param.param4
                              << " param5=" << param.param5
                              << " topDist=" << topDist[i][j] << std::endl;
                              
                }
                */
                list_vdw_.push_back(param);
            }
        }
    }

    // 计算二面角参数
    for (const auto& bond : rbond) {
        RDKit::Atom* atom1 = protac_->getAtomWithIdx(bond.first);
        RDKit::Atom* atom2 = protac_->getAtomWithIdx(bond.second);
        
        if (atom1 && atom2) {
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = protac_->getAtomNeighbors(atom1);

            while (nbrIdx != endNbrs) {
                RDKit::Atom* a = protac_->getAtomWithIdx(*nbrIdx);
                if (a->getIdx() != static_cast<unsigned int>(bond.second)) {
                    RDKit::ROMol::ADJ_ITER nbrIdx2, endNbrs2;
                    boost::tie(nbrIdx2, endNbrs2) = protac_->getAtomNeighbors(atom2);

                    while (nbrIdx2 != endNbrs2) {
                        RDKit::Atom* b = protac_->getAtomWithIdx(*nbrIdx2);
                        if (b->getIdx() != static_cast<unsigned int>(bond.first)) {
                            // 获取扭转参数
                            unsigned int torType;
                            RDKit::MMFF::MMFFTor torsionParams;
                            if (mmffProps.isValid() && mmffProps.getMMFFTorsionParams(*protac_, 
                                                                  a->getIdx(), 
                                                                  bond.first,
                                                                  bond.second, 
                                                                  b->getIdx(),
                                                                  torType,
                                                                  torsionParams)) {

                                // 打印扭转参数
                                /*
                                if(print_info){
                                    
                                    std::cout << std::fixed << std::setprecision(6)
                                              << "Dihe_params: atoms=" 
                                              << a->getIdx() << "-" << bond.first << "-" 
                                              << bond.second << "-" << b->getIdx() 
                                              << " V1=" << torsionParams.V1 
                                              << " V2=" << torsionParams.V2 
                                              << " V3=" << torsionParams.V3 << std::endl;
                                              
                                }*/

                                DihedralParam param;
                                param.atom1_idx = a->getIdx();
                                param.atom2_idx = bond.first;
                                param.atom3_idx = bond.second;
                                param.atom4_idx = b->getIdx();
                                param.v1 = torsionParams.V1;
                                param.v2 = torsionParams.V2;
                                param.v3 = torsionParams.V3;
                                //param.lowEnerDiheAngles = scanTorsion(param.v1, param.v2, param.v3);
                                //param.lowEnerDiheAngles = listTorsion(5.0);
                                DihedralEnergy lowEDihe(param.v1, param.v2, param.v3);
                                std::vector<double> minima_dihes = lowEDihe.findLocalMinima(-M_PI, M_PI, 1e-3);
                                param.lowEnerDiheAngles = lowEDihe.augmentMinima(minima_dihes, -M_PI, M_PI, static_cast<size_t>(36));
                                list_dihe_.push_back(param);
                            }
                        }
                        ++nbrIdx2;  // 移动到下一个邻居
                    }
                }
                ++nbrIdx;  // 移动到下一个邻居
            }
        }
    }
}


// 计算两个构象之间的RMSD，只考虑弹头原子
double Protac::computeRMSD(const RDKit::ROMol& mol, int confId1, int confId2) {
    // 获取两个构象
    const RDKit::Conformer& conf1 = mol.getConformer(confId1);
    const RDKit::Conformer& conf2 = mol.getConformer(confId2);
    
    // 只计算弹头原子的RMSD
    double sum_squared_dist = 0.0;
    int count = 0;
    
    for (int atom_idx : warhead_atoms_) {
        if (atom_idx < static_cast<int>(mol.getNumAtoms())) {
            const RDGeom::Point3D& pos1 = conf1.getAtomPos(atom_idx);
            const RDGeom::Point3D& pos2 = conf2.getAtomPos(atom_idx);
            
            double dx = pos1.x - pos2.x;
            double dy = pos1.y - pos2.y;
            double dz = pos1.z - pos2.z;
            
            sum_squared_dist += dx*dx + dy*dy + dz*dz;
            count++;
        }
    }
    
    if (count == 0) {
        return 0.0;  // 避免除零错误
    }
    
    return std::sqrt(sum_squared_dist / count);
}

// 聚类构象
std::vector<std::vector<int>> Protac::clusterConformers(const RDKit::ROMol& mol, double rmsdThreshold) {
    std::vector<std::vector<int>> clusters;
    unsigned int numConfs = mol.getNumConformers();
    std::vector<bool> assigned(numConfs, false);

    for (unsigned int i = 0; i < numConfs; ++i) {
        if (assigned[i]) continue;
        std::vector<int> cluster;
        cluster.push_back(i);
        assigned[i] = true;
        for (unsigned int j = i + 1; j < numConfs; ++j) {
            if (assigned[j]) continue;
            double rmsd = computeRMSD(mol, i, j);
            if (rmsd < rmsdThreshold) {
                cluster.push_back(j);
                assigned[j] = true;
            }
        }
        clusters.push_back(cluster);
    }
    return clusters;
}

void Protac::output(RDKit::SDWriter& w, 
                   std::ostream& fpro, 
                   int nKeep,
                   bool fpro_w, 
                   double rmsd_cutoff) {
    try {
        // 按能量排序
        std::sort(solutions_.begin(), solutions_.end(),
                  [](const Solution& a, const Solution& b) {
                      return a.energy < b.energy;
                  });
        
        // 创建一个临时分子来保存所有的构象
        std::unique_ptr<RDKit::RWMol> multi_conf_mol(new RDKit::RWMol(*protac_));
        std::vector<std::pair<int, size_t>> conf_to_solution; // 保存每个构象对应的solution索引
        
        // 创建所有构象
        for (size_t i = 0; i < solutions_.size(); ++i) {
            const auto& solution = solutions_[i];
            
            // 规范化二面角值到[-180, 180]范围
            std::vector<double> normalized_dihedrals;
            for (double angle : solution.dihedrals) {
                normalized_dihedrals.push_back(normalize_angle(angle));
            }
            
            // 创建新构象
            RDKit::Conformer conf(multi_conf_mol->getConformer());
            conf.setId(i);  // 设置构象ID为solution索引
            
            // 设置二面角
            bool valid_conformation = true;
            for (size_t j = 0; j < rot_dihe_.size(); ++j) {
                try {
                    const auto& dihe_j = rot_dihe_[j];
                    MolTransforms::setDihedralDeg(conf,
                                                       dihe_j[0],
                                                       dihe_j[1],
                                                       dihe_j[2],
                                                       dihe_j[3],
                                                       normalized_dihedrals[j]);
                } catch (const std::exception& e) {
                    valid_conformation = false;
                    break;
                }
            }
            
            // 验证构象中的原子距离是否合理
            if (valid_conformation) {
                multi_conf_mol->addConformer(new RDKit::Conformer(conf), true);
                conf_to_solution.push_back(std::make_pair(conf.getId(), i));
            }
        }
        
        // 聚类构象
        std::vector<std::vector<int>> clusters = clusterConformers(*multi_conf_mol, rmsd_cutoff);
        
        // 从每个集群选择能量最低的构象
        std::vector<size_t> representative_solutions;
        for (const auto& cluster : clusters) {
            if (cluster.empty()) continue;
            
            // 在当前集群中找到能量最低的构象
            int best_conf_id = cluster[0];
            double best_energy = std::numeric_limits<double>::max();
            for (int conf_id : cluster) {
                auto it = std::find_if(conf_to_solution.begin(), conf_to_solution.end(),
                                      [conf_id](const auto& p) { return p.first == conf_id; });
                if (it != conf_to_solution.end()) {
                    size_t solution_idx = it->second;
                    if (solution_idx < solutions_.size() && solutions_[solution_idx].energy < best_energy) {
                        best_energy = solutions_[solution_idx].energy;
                        best_conf_id = conf_id;
                    }
                }
            }
            
            // 添加到代表构象列表
            auto it = std::find_if(conf_to_solution.begin(), conf_to_solution.end(),
                                  [best_conf_id](const auto& p) { return p.first == best_conf_id; });
            if (it != conf_to_solution.end()) {
                representative_solutions.push_back(it->second);
            }
        }
        
        // 按能量重新排序代表构象
        std::sort(representative_solutions.begin(), representative_solutions.end(),
                 [this](size_t i, size_t j) {
                     return solutions_[i].energy < solutions_[j].energy;
                 });
        
        // 输出前nKeep个代表构象
        int output_count = 0;
        for (size_t rep_idx : representative_solutions) {
            if (output_count >= nKeep) break;
            
            const auto& solution = solutions_[rep_idx];
            
            // 创建分子的副本进行操作
            std::unique_ptr<RDKit::ROMol> protac_copy(new RDKit::ROMol(*protac_));
            RDKit::Conformer& conf = protac_copy->getConformer();
            
            // 规范化二面角值到[-180, 180]范围
            std::vector<double> normalized_dihedrals;
            for (double angle : solution.dihedrals) {
                normalized_dihedrals.push_back(normalize_angle(angle));
            }
            
            // 设置二面角
            for (size_t j = 0; j < rot_dihe_.size(); ++j) {
                const auto& dihe_j = rot_dihe_[j];
                MolTransforms::setDihedralDeg(conf,
                                                  dihe_j[0],
                                                  dihe_j[1],
                                                  dihe_j[2],
                                                  dihe_j[3],
                                                  normalized_dihedrals[j]);
            }
            
            // 添加能量分数作为分子属性
            protac_copy->setProp("SCORE", solution.energy);
            protac_copy->setProp("RANK", output_count + 1);
            
            // 写入分子
            try {
                w.write(*protac_copy);
            } catch (const std::exception& e) {
                std::cerr << "Error writing molecule: " << e.what() << std::endl;
                continue;
            }
            
            // 如果不需要输出蛋白质，直接跳过蛋白质相关计算
            if (!fpro_w) {
                output_count++;
                continue;
            }

            // 获取ref坐标（只在需要蛋白质输出时计算）
            Coords ref;
            for (int idx : idx_) {
                if (idx >= static_cast<int>(protac_copy->getNumAtoms())) {
                    std::cerr << "Invalid atom index: " << idx << std::endl;
                    continue;
                }
                const RDGeom::Point3D& pos = conf.getAtomPos(idx);
                ref.push_back({pos.x, pos.y, pos.z});
            }
            
            // 对齐蛋白质坐标
            try {
                Coords coords_pro = Align(protein_.coords, coord_subs_var, ref);
                fpro << "MODEL" << std::setw(9) << (output_count + 1) << "\n";
                fpro << "REMARK   1 SCORE " << std::fixed << std::setprecision(3) 
                     << solution.energy << "\n";
                size_t atom_index = 0;
                for (const auto& line : protein_.struct_data) {
                    if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
                        if (atom_index >= coords_pro.size()) {
                            throw std::runtime_error("Coordinate index out of bounds");
                        }
                        fpro << line.substr(0, 30)
                             << std::fixed << std::setw(8) << std::setprecision(3)
                             << coords_pro[atom_index][0]
                             << std::setw(8) << coords_pro[atom_index][1]
                             << std::setw(8) << coords_pro[atom_index][2]
                             << line.substr(54) << "\n";
                        atom_index++;
                    } else {
                        fpro << line << "\n";
                    }
                }
                fpro << "ENDMDL\n";
                
            } catch (const std::exception& e) {
                std::cerr << "Error in protein alignment: " << e.what() << std::endl;
                continue;
            }
            
            output_count++;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error in output: " << e.what() << std::endl;
        throw;
    }
}

Protac::Solution Protac::sample_single(bool verbose) {
    std::vector<double> dihe(rot_dihe_.size());
    static std::mt19937 gen(std::random_device{}());
    
    // Create a copy of the protac molecule once
    RDKit::ROMol mol_copy(*protac_);
    
    std::uniform_int_distribution<int> dist(0, 35);
    for (size_t i = 0; i < rot_dihe_.size(); ++i) {
        dihe[i] = list_dihe_[i].lowEnerDiheAngles[dist(gen)];
    }

    // Pass the copy to thread_safe_score with verbose parameter
    double ener = thread_safe_score(dihe, &mol_copy, verbose);
    return Solution{dihe, ener, std::vector<double>{}};
}

void Protac::sample(int ntotal, int nsolu, bool verbose) {
    // 创建任务vector
    if (verbose) {
        std::cout << "Ntotal: " << ntotal << std::endl;
    }
    std::vector<std::future<Solution>> futures;
    futures.reserve(ntotal);
    
    // 创建线程池
    int num_threads = std::min(processes_, static_cast<int>(std::thread::hardware_concurrency()));
    std::vector<std::thread> thread_pool;
    
    // 任务队列和结果
    std::mutex mutex;
    std::vector<Solution> solutions;
    int tasks_remaining = ntotal;
    int tasks_completed = 0;
    
    // 线程工作函数 - 移除进度条相关代码
    auto worker = [this, &mutex, &solutions, &tasks_remaining, &tasks_completed, verbose]() {
        while (true) {
            // 检查是否还有任务
            {
                std::lock_guard<std::mutex> lock(mutex);
                if (tasks_remaining <= 0) {
                    break;
                }
                --tasks_remaining;
            }
            
            // 执行采样
            Solution result = sample_single(verbose);
            
            // 保存结果
            {
                std::lock_guard<std::mutex> lock(mutex);
                solutions.push_back(result);
                ++tasks_completed;
            }
        }
    };
    
    // 启动工作线程
    for (int i = 0; i < num_threads; ++i) {
        thread_pool.emplace_back(worker);
    }
    
    // 等待所有线程完成
    for (auto& thread : thread_pool) {
        thread.join();
    }
    
    // 按能量排序
    std::sort(solutions.begin(), solutions.end(),
              [](const Solution& a, const Solution& b) {
                  return a.energy < b.energy;
              });
    
    // 选择前nsolu个解
    solutions_ = std::vector<Solution>(
        solutions.begin(),
        solutions.begin() + std::min(nsolu, static_cast<int>(solutions.size()))
    );
    
    // 执行搜索优化
    search(verbose);
}

Protac::Solution Protac::powell_minimize(const std::vector<double>& initial_guess, RDKit::ROMol* mol_copy, double tol) {
    // 规范化初始猜测值
    std::vector<double> normalized_guess = initial_guess;
    for (double& angle : normalized_guess) {
        angle = normalize_angle(angle);
    }
    
    Solution current{normalized_guess, thread_safe_score(normalized_guess, mol_copy), std::vector<double>{}};
    const int max_iter = 100;
    const int n = initial_guess.size();
    
    try {
        for (int iter = 0; iter < max_iter; ++iter) {
            double initial_energy = current.energy;
            
            for (int i = 0; i < n; ++i) {
                const double phi = (1 + std::sqrt(5.0)) / 2;
                double a = -2.0;
                double b = 2.0;
                double best_energy = current.energy;
                double best_step = 0.0;

                for (int j = 0; j < 10; ++j) {
                    double x1 = normalize_angle(b - (b - a) / phi);
                    double x2 = normalize_angle(a + (b - a) / phi);

                    std::vector<double> guess1 = current.parameters;
                    std::vector<double> guess2 = current.parameters;
                    
                    if (i < static_cast<int>(guess1.size()) && i < static_cast<int>(guess2.size())) {
                        guess1[i] = normalize_angle(current.parameters[i] + x1);
                        guess2[i] = normalize_angle(current.parameters[i] + x2);
                        
                        for (double& angle : guess1) angle = normalize_angle(angle);
                        for (double& angle : guess2) angle = normalize_angle(angle);
                        
                        double energy1 = thread_safe_score(guess1, mol_copy);
                        double energy2 = thread_safe_score(guess2, mol_copy);

                        if (energy1 < energy2) {
                            b = normalize_angle(x2);
                            if (energy1 < best_energy) {
                                best_energy = energy1;
                                best_step = normalize_angle(x1);
                            }
                        } else {
                            a = normalize_angle(x1);
                            if (energy2 < best_energy) {
                                best_energy = energy2;
                                best_step = normalize_angle(x2);
                            }
                        }
                    }
                }

                if (i < static_cast<int>(current.parameters.size())) {
                    current.parameters[i] = normalize_angle(current.parameters[i] + best_step);
                    for (double& angle : current.parameters) {
                        angle = normalize_angle(angle);
                    }
                    current.energy = best_energy;
                }
            }
            
            if (std::abs(initial_energy - current.energy) < tol) {
                break;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error in powell_minimize: " << e.what() << std::endl;
    }
    
    for (double& angle : current.parameters) {
        angle = normalize_angle(angle);
    }
    
    return current;
}

Protac::Solution Protac::search_single(const Solution& initial_solution) {
    // Create a copy of the protac molecule once for this search
    RDKit::ROMol mol_copy(*protac_);
    
    // 第一阶段参数
    std::vector<double> temperatures = {60.0, 50.0, 40.0, 30.0, 20.0};
    double alpha = 1.1;
    // 保存最佳解
    Solution best = initial_solution;    
    // 第一阶段：高温搜索
    for (double T : temperatures) {
        alpha *= 0.9;
        // Powell最小化
        Solution minimized = powell_minimize(best.dihedrals, &mol_copy, 0.01);
        if (minimized.energy < best.energy) {
            best = minimized;
        }
        
        Solution prev = best;
        // Monte Carlo搜索
        for (int i = 0; i < 200; ++i) {
            // 生成随机扰动
            std::vector<double> delta(rot_dihe_.size());
            std::uniform_real_distribution<double> perturbation(-10.0, 10.0);
            
            for (size_t j = 0; j < rot_dihe_.size(); ++j) {
                delta[j] = perturbation(rng_) * alpha;
            }
            
            // 计算新构象
            std::vector<double> current = prev.dihedrals;
            for (size_t j = 0; j < current.size(); ++j) {
                current[j] += delta[j];
            }
            
            // 计算能量
            double score = this->thread_safe_score(current, &mol_copy);
            // Metropolis准则
            if (score < prev.energy) {
                prev = Solution{current, score, std::vector<double>{}};
                if (score < best.energy) {
                    best = prev;
                }
            } else {
                double probability = std::exp((prev.energy - score) / T);
                if (uniform_dist_(rng_) < probability) {
                    prev = Solution{current, score, std::vector<double>{}};
                }
            }
        }
    }
    
    // 第二阶段：低温精细搜索
    {
        Solution minimized = powell_minimize(best.dihedrals, &mol_copy, 0.01);
        if (minimized.energy < best.energy) {
            best = minimized;
        }
        
        Solution prev = best;
        const double T = 10.0;
        
        for (int i = 0; i < 300; ++i) {
            // 生成小幅随机扰动
            std::vector<double> delta(rot_dihe_.size());
            std::uniform_real_distribution<double> perturbation(-5.0, 5.0);
            
            for (size_t j = 0; j < rot_dihe_.size(); ++j) {
                delta[j] = perturbation(rng_);
            }
            
            // 计算新构象
            std::vector<double> current = prev.dihedrals;
            for (size_t j = 0; j < current.size(); ++j) {
                current[j] += delta[j];
            }
            
            // 计算能量
            double score = this->thread_safe_score(current, &mol_copy);
            // Metropolis准则
            if (score < prev.energy) {
                prev = Solution{current, score, std::vector<double>{}};
                if (score < best.energy) {
                    best = prev;
                }
            } else {
                double probability = std::exp((prev.energy - score) / T);
                if (uniform_dist_(rng_) < probability) {
                    prev = Solution{current, score, std::vector<double>{}};
                }
            }
        }
    }
    
    // 最终优化
    Solution final = powell_minimize(best.dihedrals, &mol_copy, 0.01);
    if (final.energy < best.energy) {
        best = final;
    }
    return best;
}

ThreadPool::ThreadPool(size_t num_threads) : stop(false) {
    for (size_t i = 0; i < num_threads; ++i) {
        workers.emplace_back([this] {
            for (;;) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                    if (this->stop && this->tasks.empty()) return;
                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                }
                task();
            }
        });
    }
}

ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for (std::thread &worker : workers) {
        worker.join();
    }
}

template<class F>
auto ThreadPool::enqueue(F&& f) -> std::future<typename std::result_of<F()>::type> {
    using return_type = typename std::result_of<F()>::type;

    auto task = std::make_shared<std::packaged_task<return_type()>>(std::forward<F>(f));
    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        tasks.emplace([task]() { (*task)(); });
    }
    condition.notify_one();
    return res;
}

// 在 Protac::search() 中使用 ThreadPool
void Protac::search(bool verbose) {
    int num_threads = std::min(processes_, static_cast<int>(std::thread::hardware_concurrency()));
    ThreadPool pool(num_threads); // 使用 num_threads 创建线程池
    std::vector<std::future<Solution>> futures;
    futures.reserve(solutions_.size());
    
    // 提交所有任务
    for (const auto& solution : solutions_) {
        futures.push_back(pool.enqueue(
            [this, solution]() { return search_single(solution); }
        ));
    }

    // 准备进度条（仅在verbose模式下）
    std::unique_ptr<ProgressBar> progress_bar;
    if (verbose) {
        progress_bar = std::make_unique<ProgressBar>(futures.size());
    }
    
    // 收集结果
    solutions_.clear();
    for (size_t i = 0; i < futures.size(); ++i) {
        solutions_.push_back(futures[i].get());
        if (verbose && progress_bar) {
            progress_bar->update(i + 1); // 更新进度条，传递当前已完成的任务数
        }
    }

    if (verbose && progress_bar) {
        progress_bar->finish(); // 结束进度条
    }
}

double Protac::e_intra(const RDKit::ROMol* mol) const {
    const RDKit::Conformer& conf = mol->getConformer();    
    
    double vdw = 0.0;
    
    std::vector<RDGeom::Point3D> coors;
    coors.reserve(mol->getNumAtoms());
    for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
        coors.push_back(conf.getAtomPos(i));
    }
    for (const auto& p : list_vdw_) {
        double dist = (coors[p.atom1_idx] - coors[p.atom2_idx]).length();
        double term1 = p.param2 / (dist + p.param3);
        double term1_pow = std::pow(term1, 7);
        double term2 = p.param4 / (std::pow(dist, 7) + p.param5);
        double result = term1_pow * (term2 - 2.0);
        //std::cout << "[VDW]: " << p.param1 * result << ": " << p.atom1_idx << "-" << p.atom2_idx << " " << dist << std::endl;
        vdw += (p.param1 * result > 3.0) ? 3.0 : p.param1 * result;
    }
    double dihe = 0.0;
    for (const auto& d : list_dihe_) {
        double angle = MolTransforms::getDihedralRad(conf,
                                                   d.atom1_idx,
                                                   d.atom2_idx,
                                                   d.atom3_idx,
                                                   d.atom4_idx);
        
        dihe += 0.5 * (
            d.v1 * (1.0 + std::cos(angle)) +
            d.v2 * (1.0 - std::cos(2.0 * angle)) +
            d.v3 * (1.0 + std::cos(3.0 * angle))
        );
    }

    return vdw + dihe;
}



double Protac::thread_safe_score(const std::vector<double>& dihe, RDKit::ROMol* mol_copy, bool print_info) {
    // Use the provided molecule copy instead of creating a new one
    RDKit::Conformer& conf = mol_copy->getConformer();
    
    // 首先规范化所有二面角
    std::vector<double> normalized_dihe = dihe;
    for (double& angle : normalized_dihe) {
        angle = normalize_angle(angle);
    }
    
    // 对每个二面角
    for (size_t j = 0; j < rot_dihe_.size(); ++j) {
        const auto& atoms = rot_dihe_[j];
        // 设置二面角
        MolTransforms::setDihedralDeg(conf,
                                    atoms[0],
                                    atoms[1],
                                    atoms[2],
                                    atoms[3],
                                    normalized_dihe[j]);
    }

    // 获取所有protac原子位置
    std::vector<RDGeom::Point3D> positions;
    positions.reserve(mol_copy->getNumAtoms());
    for (size_t i = 0; i < mol_copy->getNumAtoms(); ++i) {
        positions.push_back(conf.getAtomPos(i));
    }

    // 计算分子内能量
    double e_in = e_intra(mol_copy) - E_intra_ref_;
    if (e_in > paras["ub_strain"]) {
        e_in = paras["ub_strain"];
    }

    // 计算q_flex与锚定蛋白的相互作用能
    double e_anchor = 0.0;
    Coords coors_flex;
    for (const auto& q : q_flex_) {
        coors_flex.push_back({
            positions[std::get<0>(q).value()].x,
            positions[std::get<0>(q).value()].y,
            positions[std::get<0>(q).value()].z
        });
    }
    e_anchor = GetGridEn(grid_anchor_, coors_flex, q_flex_);
    
    // 蛋白质-蛋白质相互作用
    double e_pp = 0.0;
    Coords ref;
    for (int idx : idx_) {
        const auto& pos = positions[idx];
        ref.push_back({pos.x, pos.y, pos.z});
    }
    
    Coords coords_pro = Align(protein_.coords, coord_subs_var, ref);
    e_pp = GetGridEn(grid_anchor_, coords_pro, protein_.para);

    // q_anchor与柔性蛋白的相互作用能
    double e_flex = 0.0;
    Coords coors_anchor;
    for (const auto& q : q_anchor_) {
        coors_anchor.push_back({
            positions[std::get<0>(q).value()].x,
            positions[std::get<0>(q).value()].y,
            positions[std::get<0>(q).value()].z
        });
    }
        
    Coords aligned_coors = Align2(coors_anchor, coord_subs_var, ref, translation_);
    e_flex = GetGridEn(grid_flex_, aligned_coors, q_anchor_);

    return e_in + e_anchor + e_pp + e_flex;
}

void Protac::printProtacInfo(){
    std::cout << "Number of atoms: " << protac_->getNumAtoms() << std::endl;
    std::cout << "Number of bonds: " << protac_->getNumBonds() << std::endl;
    std::cout << "IntraEnergy    : " << e_intra(protac_.get()) << std::endl;
    // 将MolBlock写入到文件
    std::ofstream outFile("initialized_protac.sdf");
    if (outFile.is_open()) {
        outFile << RDKit::MolToMolBlock(*protac_) << "\n$$$$\n";
        outFile.close();
    } else {
        std::cerr << "Error: Unable to open file 'init_protac.sdf' for writing." << std::endl;
    }
    for (size_t i = 0; i < rot_dihe_.size(); ++i) {
        const auto& dihedral = rot_dihe_[i];
        RDKit::Conformer& conf = protac_->getConformer();
        
        // 获取原子信息
        RDKit::Atom* a1 = protac_->getAtomWithIdx(dihedral[0]);
        RDKit::Atom* a2 = protac_->getAtomWithIdx(dihedral[1]);
        RDKit::Atom* a3 = protac_->getAtomWithIdx(dihedral[2]);
        RDKit::Atom* a4 = protac_->getAtomWithIdx(dihedral[3]);
        
        // 获取原子位置
        const RDGeom::Point3D& p1 = conf.getAtomPos(dihedral[0]);
        const RDGeom::Point3D& p2 = conf.getAtomPos(dihedral[1]);
        const RDGeom::Point3D& p3 = conf.getAtomPos(dihedral[2]);
        const RDGeom::Point3D& p4 = conf.getAtomPos(dihedral[3]);
        
        // 计算键长
        double d12 = (p1 - p2).length();
        double d23 = (p2 - p3).length();
        double d34 = (p3 - p4).length();
        
        std::cout << "Dihedral " << i << ":" << std::endl;
        std::cout << "  Atoms: " << dihedral[0] << "(" << a1->getSymbol() << ")-"
                  << dihedral[1] << "(" << a2->getSymbol() << ")-"
                  << dihedral[2] << "(" << a3->getSymbol() << ")-"
                  << dihedral[3] << "(" << a4->getSymbol() << ")" << std::endl;
        std::cout << "  Bond lengths: " << d12 << ", " << d23 << ", " << d34 << " Å" << std::endl;
        std::cout << "  Positions:" << std::endl;
        std::cout << "    " << dihedral[0] << ": (" << p1.x << ", " << p1.y << ", " << p1.z << ")" << std::endl;
        std::cout << "    " << dihedral[1] << ": (" << p2.x << ", " << p2.y << ", " << p2.z << ")" << std::endl;
        std::cout << "    " << dihedral[2] << ": (" << p3.x << ", " << p3.y << ", " << p3.z << ")" << std::endl;
        std::cout << "    " << dihedral[3] << ": (" << p4.x << ", " << p4.y << ", " << p4.z << ")" << std::endl;
        std::cout << " END PRINT PROTAC INFO"<< std::endl;
    }
}

// 新增score_only功能实现
double Protac::score_only(bool verbose) {
    // 对当前分子构象进行评分，不改变二面角
    RDKit::ROMol mol_copy(*protac_);
    
    // 获取当前二面角值
    const RDKit::Conformer& conf = protac_->getConformer();
    std::vector<double> current_dihedrals;
    current_dihedrals.reserve(rot_dihe_.size());
    for (size_t i = 0; i < rot_dihe_.size(); ++i) {
        const auto& atoms = rot_dihe_[i];
        double angle = MolTransforms::getDihedralDeg(conf, atoms[0], atoms[1], atoms[2], atoms[3]);
        current_dihedrals.push_back(angle);
    }
    
    double total_energy = thread_safe_score(current_dihedrals, &mol_copy, true);
    
    if (verbose) {
        std::cout << "\n=== Score Only Mode ===" << std::endl;
        std::cout << "Evaluating current PROTAC conformation..." << std::endl;
        
        std::cout << "Current dihedral angles:" << std::endl;
        for (size_t i = 0; i < rot_dihe_.size(); ++i) {
            const auto& atoms = rot_dihe_[i];
            std::cout << "  Dihedral " << i << " [" << atoms[0] << "-" << atoms[1] 
                      << "-" << atoms[2] << "-" << atoms[3] << "]: " 
                      << std::fixed << std::setprecision(2) << current_dihedrals[i] << "°" << std::endl;
        }
        
        // 计算各能量组分
        std::cout << "\nEnergy components:" << std::endl;
        
        // 分子内能量
        double e_in = e_intra(protac_.get()) - E_intra_ref_;
        if (e_in > paras["ub_strain"]) {
            e_in = paras["ub_strain"];
        }
        std::cout << "  Intramolecular energy (strain): " << std::fixed << std::setprecision(3) << e_in << " kcal/mol" << std::endl;
        
        // 获取所有protac原子位置
        std::vector<RDGeom::Point3D> positions;
        positions.reserve(protac_->getNumAtoms());
        for (size_t i = 0; i < protac_->getNumAtoms(); ++i) {
            positions.push_back(conf.getAtomPos(i));
        }
        
        // q_flex与锚定蛋白的相互作用能
        Coords coors_flex;
        for (const auto& q : q_flex_) {
            coors_flex.push_back({
                positions[std::get<0>(q).value()].x,
                positions[std::get<0>(q).value()].y,
                positions[std::get<0>(q).value()].z
            });
        }
        double e_anchor = GetGridEn(grid_anchor_, coors_flex, q_flex_);
        std::cout << "  Flexible warhead - Anchor protein: " << std::fixed << std::setprecision(3) << e_anchor << " kcal/mol" << std::endl;
        
        // 蛋白质-蛋白质相互作用
        Coords ref;
        for (int idx : idx_) {
            const auto& pos = positions[idx];
            ref.push_back({pos.x, pos.y, pos.z});
        }
        Coords coords_pro = Align(protein_.coords, coord_subs_var, ref);
        double e_pp = GetGridEn(grid_anchor_, coords_pro, protein_.para);
        std::cout << "  Protein-protein interaction: " << std::fixed << std::setprecision(3) << e_pp << " kcal/mol" << std::endl;
        
        // q_anchor与柔性蛋白的相互作用能
        Coords coors_anchor;
        for (const auto& q : q_anchor_) {
            coors_anchor.push_back({
                positions[std::get<0>(q).value()].x,
                positions[std::get<0>(q).value()].y,
                positions[std::get<0>(q).value()].z
            });
        }
        Coords aligned_coors = Align2(coors_anchor, coord_subs_var, ref, translation_);
        double e_flex = GetGridEn(grid_flex_, aligned_coors, q_anchor_);
        std::cout << "  Anchor warhead - Flexible protein: " << std::fixed << std::setprecision(3) << e_flex << " kcal/mol" << std::endl;
        
        std::cout << "  ------------------------" << std::endl;
        std::cout << "  Total energy: " << std::fixed << std::setprecision(3) << total_energy << " kcal/mol" << std::endl;
        std::cout << "=======================" << std::endl;
    }
    
    return total_energy;
}

double Protac::score_only(const std::vector<double>& dihe, bool verbose) {
    if (dihe.size() != rot_dihe_.size()) {
        throw std::runtime_error("Dihedral vector size (" + std::to_string(dihe.size()) + 
                               ") does not match number of rotatable dihedrals (" + 
                               std::to_string(rot_dihe_.size()) + ")");
    }
    
    // 创建分子副本并设置指定的二面角
    RDKit::ROMol mol_copy(*protac_);
    RDKit::Conformer& conf = mol_copy.getConformer();
    
    // 设置用户指定的二面角
    for (size_t j = 0; j < rot_dihe_.size(); ++j) {
        const auto& atoms = rot_dihe_[j];
        MolTransforms::setDihedralDeg(conf, atoms[0], atoms[1], atoms[2], atoms[3], dihe[j]);
    }
    
    double total_energy = thread_safe_score(dihe, &mol_copy, true);
    
    if (verbose) {
        std::cout << "\n=== Score Only Mode (Custom Dihedrals) ===" << std::endl;
        std::cout << "Evaluating PROTAC with specified dihedral angles..." << std::endl;
        
        std::cout << "Specified dihedral angles:" << std::endl;
        for (size_t i = 0; i < rot_dihe_.size(); ++i) {
            const auto& atoms = rot_dihe_[i];
            std::cout << "  Dihedral " << i << " [" << atoms[0] << "-" << atoms[1] 
                      << "-" << atoms[2] << "-" << atoms[3] << "]: " 
                      << std::fixed << std::setprecision(2) << dihe[i] << "°" << std::endl;
        }
        
        // 计算各能量组分（重新计算基于新构象）
        std::cout << "\nEnergy components:" << std::endl;
        
        // 获取所有protac原子位置
        std::vector<RDGeom::Point3D> positions;
        positions.reserve(mol_copy.getNumAtoms());
        for (size_t i = 0; i < mol_copy.getNumAtoms(); ++i) {
            positions.push_back(conf.getAtomPos(i));
        }
        
        // 分子内能量
        double e_in = e_intra(&mol_copy) - E_intra_ref_;
        if (e_in > paras["ub_strain"]) {
            e_in = paras["ub_strain"];
        }
        std::cout << "  Intramolecular energy (strain): " << std::fixed << std::setprecision(3) << e_in << " kcal/mol" << std::endl;
        
        // q_flex与锚定蛋白的相互作用能
        Coords coors_flex;
        for (const auto& q : q_flex_) {
            coors_flex.push_back({
                positions[std::get<0>(q).value()].x,
                positions[std::get<0>(q).value()].y,
                positions[std::get<0>(q).value()].z
            });
        }
        double e_anchor = GetGridEn(grid_anchor_, coors_flex, q_flex_);
        std::cout << "  Flexible warhead - Anchor protein: " << std::fixed << std::setprecision(3) << e_anchor << " kcal/mol" << std::endl;
        
        // 蛋白质-蛋白质相互作用
        Coords ref;
        for (int idx : idx_) {
            const auto& pos = positions[idx];
            ref.push_back({pos.x, pos.y, pos.z});
        }
        Coords coords_pro = Align(protein_.coords, coord_subs_var, ref);
        double e_pp = GetGridEn(grid_anchor_, coords_pro, protein_.para);
        std::cout << "  Protein-protein interaction: " << std::fixed << std::setprecision(3) << e_pp << " kcal/mol" << std::endl;
        
        // q_anchor与柔性蛋白的相互作用能
        Coords coors_anchor;
        for (const auto& q : q_anchor_) {
            coors_anchor.push_back({
                positions[std::get<0>(q).value()].x,
                positions[std::get<0>(q).value()].y,
                positions[std::get<0>(q).value()].z
            });
        }
        Coords aligned_coors = Align2(coors_anchor, coord_subs_var, ref, translation_);
        double e_flex = GetGridEn(grid_flex_, aligned_coors, q_anchor_);
        std::cout << "  Anchor warhead - Flexible protein: " << std::fixed << std::setprecision(3) << e_flex << " kcal/mol" << std::endl;
        
        std::cout << "  ------------------------" << std::endl;
        std::cout << "  Total energy: " << std::fixed << std::setprecision(3) << total_energy << " kcal/mol" << std::endl;
        std::cout << "=======================" << std::endl;
    }
    
    return total_energy;
}




