#include <iostream>
#include <map>
#include <deque>
#include <string>
#include "DetectorAnalyzer.h"

void PreAnalyzer(bool debugMode = true) {
    // 出力ファイル名
    std::string outputFileName = "analyzed_output.root";
    std::cout << "[INFO] Analyzer initialized. Output: " << outputFileName << std::endl;

    if (debugMode) {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "  DEBUG MODE ON: Processing first 10 minutes only." << std::endl;
        std::cout << "==========================================\n" << std::endl;
    } else {
        std::cout << "FULL MODE: Processing ALL data." << std::endl;
    }

    // 解析クラスの初期化 (タイムウィンドウ 100ns)
    DetectorAnalyzer analyzer(100, outputFileName);

    // 時間制限設定
    if (debugMode) {
        analyzer.setTimeLimitMinutes(60.0);
    } else {
        analyzer.setTimeLimitMinutes(-1.0); // 無制限
    }
analyzer.loadOfflineStripList("offline_strip_list.csv");
    // CSV読み込み (X1, Y1 -> Mod 1 / X2, Y2 -> Mod 0)
    // DetectorAnalyzer側でYオフセット(+64)が自動適用されます
    std::cout << "[INFO] Loading CSV mapping files..." << std::endl;
    try {
        analyzer.loadConversionFactorsFromCSV("X1.csv", "X1", 1);
        analyzer.loadConversionFactorsFromCSV("Y1.csv", "Y1", 1);
        analyzer.loadConversionFactorsFromCSV("X2.csv", "X2", 0);
        analyzer.loadConversionFactorsFromCSV("Y2.csv", "Y2", 0);
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Failed to load CSV files: " << e.what() << std::endl;
        return;
    }

    // ファイル検索とソート
    // ディレクトリ内の .edb ファイルを探し、モジュールごとに分類して時系列ソートします
    std::map<int, std::deque<std::string>> fileQueues;
    std::string dataDir = "./"; 

    std::cout << "[INFO] Searching and sorting binary files in " << dataDir << "..." << std::endl;
    analyzer.loadAndSortFiles(dataDir, fileQueues);

    // 解析実行
    if (!fileQueues.empty()) {
        analyzer.processBinaryFiles(fileQueues);
    } else {
        std::cerr << "[WARN] No valid data files found to process." << std::endl;
    }

    std::cout << "[INFO] PreAnalyzer completed." << std::endl;
}