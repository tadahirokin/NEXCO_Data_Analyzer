#include "DetectorAnalyzer.h"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <deque>
#include <algorithm>
#include <regex>
#include <TSystem.h> 

// ファイル名からモジュールIDを抽出
// ルール: "後ろから2番目の3桁の数値"
// 例: PSD_00_001_000.edb -> 1
int extractModuleID(const std::string& path) {
    try {
        // 正規表現: アンダースコア + 3桁数字 + アンダースコア + 数字(ファイル番号) + .edb (行末)
        std::regex re("_([0-9]{3})_[0-9]+\\.edb$");
        std::smatch match;
        if (std::regex_search(path, match, re)) {
            return std::stoi(match.str(1));
        }
    } catch (...) { }
    return -1; // マッチしない場合
}

long extractDateFromPath(const std::string& path) {
    try {
        std::regex re("PSD[^_]*_([0-9]{8})");
        std::smatch match;
        if (std::regex_search(path, match, re)) return std::stol(match.str(1));
    } catch (...) { }
    return 0; 
}

std::vector<std::string> findEdbFilesFromFolder(const std::string& folderPath) {
    std::cout << "[INFO] Searching for .edb files in: " << folderPath << std::endl;
    std::vector<std::string> fileNames;
    TSystem* system = gSystem;
    std::vector<std::string> foldersToSearch;
    foldersToSearch.push_back(folderPath);

    void* dir = system->OpenDirectory(folderPath.c_str());
    if (dir) {
        const char* entry;
        while ((entry = system->GetDirEntry(dir))) {
            std::string sEntry = entry;
            if (sEntry == "." || sEntry == "..") continue;
            std::string fullPath = folderPath + "/" + sEntry;
            FileStat_t st;
            if (system->GetPathInfo(fullPath.c_str(), st) == 0 && R_ISDIR(st.fMode)) {
                foldersToSearch.push_back(fullPath);
            }
        }
        system->FreeDirectory(dir);
    }

    for (const auto& fPath : foldersToSearch) {
        void* subDir = system->OpenDirectory(fPath.c_str());
        if (!subDir) continue;
        const char* subEntry;
        while ((subEntry = system->GetDirEntry(subDir))) {
            std::string sSubEntry = subEntry;
            if (sSubEntry.length() >= 4 && sSubEntry.substr(sSubEntry.length() - 4) == ".edb") {
                if (sSubEntry.find("._") == 0) continue;
                fileNames.push_back(fPath + "/" + sSubEntry);
            }
        }
        system->FreeDirectory(subDir);
    }
    return fileNames;
}

// debugMode = true: 10分解析
void PreAnalyzer(bool debugMode = true) {
    const int TIME_WINDOW_NS = 100;
    std::string rootFileName = "analyzed_output.root";
    std::string dataFolder = "./"; 

    DetectorAnalyzer analyzer(TIME_WINDOW_NS, rootFileName);

    if (debugMode) {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "  DEBUG MODE ON: Processing first 10 minutes only." << std::endl;
        std::cout << "  (To run FULL analysis, execute: PreAnalyzer(false))" << std::endl;
        std::cout << "==========================================\n" << std::endl;
        analyzer.setTimeLimitMinutes(10.0);
    } else {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "  FULL MODE: Processing ALL data." << std::endl;
        std::cout << "==========================================\n" << std::endl;
        analyzer.setTimeLimitMinutes(-1); 
    }

    // =========================================================
    // マッピング設定
    // =========================================================
    try {
        std::cout << "[INFO] Loading CSV mapping files..." << std::endl;
        // Mod 1 -> X1, Y1 / Mod 0 -> X2, Y2
        analyzer.loadConversionFactorsFromCSV("X1.csv", "X1", 1);
        analyzer.loadConversionFactorsFromCSV("Y1.csv", "Y1", 1);
        analyzer.loadConversionFactorsFromCSV("X2.csv", "X2", 0);
        analyzer.loadConversionFactorsFromCSV("Y2.csv", "Y2", 0);
    } catch (const std::runtime_error& e) {
        std::cerr << "[ERROR] CSV loading failed: " << e.what() << std::endl;
        return;
    }

    // ファイル検索
    std::vector<std::string> allFiles = findEdbFilesFromFolder(dataFolder);
    if (allFiles.empty()) {
        std::cerr << "[CRITICAL] No data files found." << std::endl;
        return;
    }

    // 【重要】ファイルをモジュールIDごとに分類 (Map作成)
    std::map<int, std::deque<std::string>> fileQueues;
    int unknownCount = 0;
    for (const auto& f : allFiles) {
        int modID = extractModuleID(f);
        if (modID != -1) {
            fileQueues[modID].push_back(f);
        } else {
            unknownCount++;
        }
    }

    // 分類結果の表示 & ソート
    std::cout << "------------------------------------------" << std::endl;
    std::cout << " File Classification Result:" << std::endl;
    if (fileQueues.empty()) {
        std::cerr << "[ERROR] No valid module IDs found in filenames!" << std::endl;
        if (unknownCount > 0) std::cerr << " -> " << unknownCount << " files failed to match pattern '_XXX_XXX.edb'." << std::endl;
        return;
    }
    for (auto& [mod, queue] : fileQueues) {
        // 名前順(≒時刻順)にソート
        std::sort(queue.begin(), queue.end());
        std::cout << "  Module " << mod << ": " << queue.size() << " files found." << std::endl;
        if (!queue.empty()) std::cout << "    (First: " << queue.front() << ")" << std::endl;
    }
    if (unknownCount > 0) {
        std::cout << "  (Ignored " << unknownCount << " files with unknown format)" << std::endl;
    }
    std::cout << "------------------------------------------" << std::endl;

    // 解析実行 (修正: vectorではなくmapを渡す)
    analyzer.processBinaryFiles(fileQueues);

    if (analyzer.getRawTotData().empty()) {
        std::cout << "[WARNING] No valid events found." << std::endl;
    } else {
        std::cout << "[INFO] PreAnalyzer completed. Check " << rootFileName << std::endl;
    }
}
