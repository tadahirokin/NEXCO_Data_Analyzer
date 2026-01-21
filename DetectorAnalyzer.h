#ifndef DETECTOR_ANALYZER_H
#define DETECTOR_ANALYZER_H

#include <string>
#include <vector>
#include <map>
#include <deque>
#include <set>
#include <utility>
#include <fstream>
#include <iostream>
#include <climits> 
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMultiGraph.h> 
#include <TCanvas.h>
#include <TArrow.h> 
#include <TLine.h> 
#include <TBox.h>
#include <TSpectrum.h>
#include <TROOT.h>
#include <TStyle.h>
#include <chrono>
#include <iomanip>

// ID定義
const int ID_X1 = 0;
const int ID_Y1 = 1;
const int ID_X2 = 2;
const int ID_Y2 = 3;

// システム定数
const int MAX_MODULES = 16;   
const int MAX_SYS_CH = 512;
const int Y_CHANNEL_OFFSET = 64; 

// 閾値設定
const int TOT_THRESHOLD_NS = 1000; 



const unsigned long long ANALYSIS_START_SKIP_NS = 180000000000ULL; // 3分 (3 * 60 * 10^9)
// 解析ウィンドウ設定 (10分)
const unsigned long long GAIN_ANALYSIS_WINDOW_NS = 600000000000ULL; 
const unsigned long long CHANNEL_TIMEOUT_NS = 5000000000ULL;       // 5秒

// ゲインモニタリング解像度 (20ns/bin)
const int MONITOR_HIST_BINS = 4000;
const double MONITOR_HIST_MAX_TOT = 100000.0;
const double BIN_WIDTH_NS = 20.0; 

// ピーク探索用初期範囲
const double PEAK_SEARCH_MIN_TOT = 30000.0;
const double PEAK_SEARCH_MAX_TOT = 80000.0;
const int SILENCE_THRESHOLD_BINS = 5; // 連続ゼロ判定

// 互換性のために必要
const int COVELL_WINDOW_NS = 1000;

// マッチング探索範囲
const int MATCHING_SEARCH_RANGE_BINS = 50; 
const int SHIFT_CALC_ERROR = -999999;

class DetectorAnalyzer {
public:
    struct Event {
        int type;
        int strip;
        int module;
        unsigned long long eventTime_ns;
        int tot;
        int sysCh;
    };

    struct ChConfig {
        int detTypeID;
        int strip; 
        std::string cachedName;
        int pIndex;
    };

    struct LUTItem {
        bool isValid;
        int detTypeID;
        int strip;
        int pIndex;
    };

    DetectorAnalyzer(int timeWindow_ns, const std::string& outputFileName);
    ~DetectorAnalyzer();

    // ★モード切替スイッチ
    void setMuonAnalysisMode(bool enable);

    void loadConversionFactorsFromCSV(const std::string& csvFileName, const std::string& detName, int moduleID);
    void loadOfflineStripList(const std::string& csvFileName);
    void loadAndSortFiles(const std::string& directory, std::map<int, std::deque<std::string>>& fileQueues);
    void setTimeLimitMinutes(double minutes);
    void processBinaryFiles(std::map<int, std::deque<std::string>>& fileQueues);

    void printSearchStatus();
    bool isAnalysisStarted_; 

private:
    void setupTree();
    void setupIdealTree(); // ミュオン解析用Treeの初期化ヘルパー

    // ★モード管理フラグ
    bool enableMuonAnalysis_;

    // ★処理の実装分離
    bool processChunk_Fast(const std::vector<Event>& sortedEvents); // ゲイン解析のみ (爆速)
    bool processChunk_Muon(const std::vector<Event>& sortedEvents); // 自動抽出あり (高負荷)
    bool processChunk(const std::vector<Event>& sortedEvents);      // エントリポイント

    // 理想イベント抽出用ファイル・Tree
    TFile* idealFile_;
    TTree* idealTree_;

    // idealTree用バッファ
    int    b_i_runID;
    int    b_i_eventID;
    int    b_i_nHits;
    double b_i_width;
    std::vector<int>    *b_i_type;
    std::vector<int>    *b_i_strip;
    std::vector<int>    *b_i_tot;
    std::vector<double> *b_i_time;

    // 抽出アルゴリズム用バッファ
    std::deque<Event> analysisBuffer_;
    int globalEventID_Ideal_;

    // 幅ヒストグラム
    TH1F *hEventWidth_All;
    TH1F *hEventWidth_CondA;
    TH1F *hEventWidth_CondB;

    // アルゴリズム用メソッド
    void processEventExtraction(); 
    bool checkConditionA(const std::deque<Event>& buf, size_t startIdx, double& outWidth);
    bool checkConditionB(const std::deque<Event>& buf, size_t startIdx, double& outWidth);
    bool hasAdjacentPair(const std::vector<int>& strips);
    bool isAdjacentToOffline(int strip, int detTypeID);

    // PDF生成
    void generatePDFReport();
    
    // 同期・読み込み
    bool syncDataStream(std::ifstream& ifs, long long& skippedBytes);
    std::pair<unsigned long long, bool> readEventsFromFile(const std::string& fileName, int modID, std::vector<Event>& rawEvents, long long& offset);
    
    unsigned long long getSafeTime(const std::map<int, unsigned long long>& lastTimes);
    std::string getRunSignature(const std::string& fileName);
    short getRunID(const std::string& prefix);
    std::string parseRunPrefix(const std::string& fullPath);
    
    // ゲイン解析ロジック
    bool analyzeGainShift(unsigned long long currentTime_ns);
    int findRightMostPeak(TH1F* hist); 
    void determineIntegrationRange(TH1F* hist, int peakBin, int& outMinBin, int& outMaxBin);
    
    double calculatePeakAreaCovell(TH1F* hist, double peakPos);
    double calculateResidual(TH1F* hTarget, TH1F* hTemplate, int shiftBins, int binMin, int binMax);
    int findBestShift(TH1F* hTarget, TH1F* hTemplate, int binMin, int binMax, const std::string& debugLabel);
    
    void printFileTail(const std::string& filePath, long long currentOffset);
    void calculateEffectiveTime();
    void printLog(const std::string& msg) { std::cout << msg << std::endl; }

    // --- メンバ変数 ---
    int timeWindow_ns_;
    TFile* outputFile_;
    TTree* tree_;
    TH1F *hDeltaT_Nearest, *hDeltaT_n8, *hGlobalPIndexMap;
    TGraph* gT0Check;
    int t0GlobalCount_;
    unsigned long long firstT0Time_;

    // TTreeバッファ
    Char_t   b_type;   
    Char_t   b_strip;  
    Int_t    b_tot;    
    Long64_t b_time;   

    // ゲイン追従・解析用
    std::map<std::pair<int, int>, double> initialPeakPosMap_;    
    std::map<std::pair<int, int>, double> cumulativeShiftNsMap_; 
    std::map<std::pair<int, int>, int>    rangeMinMap_;          
    std::map<std::pair<int, int>, int>    rangeMaxMap_;          

    unsigned long long analysisDuration_ns_;
    unsigned long long globalStartTimestamp_;
    unsigned long long currentEventTime_ns_;
    unsigned long long lastGlobalAnalysisTime_ns_; 
    
    bool isFirstGainAnalysis_;
    bool isCsvHeaderWritten_;
    bool isTemplateCaptured_;

    unsigned long long totalDataSize_;
    unsigned long long processedDataSize_;
    std::chrono::system_clock::time_point analysisStartTime_;
    bool dashboardShown_;

    std::map<std::string, std::map<int, std::pair<int, int>>> nameStripToMapKey_;
    std::map<std::pair<int, int>, std::set<int>> foundChannels_; 
    
    std::map<int, std::string> currentRunPrefix_; 
    std::set<std::pair<int, int>> ignoredMonitorChannels_; 
    std::set<std::pair<int, int>> offlineStrips_;

    std::map<int, unsigned long long> channelLastRefHitTime_; 
    unsigned long long lastTimeoutCheckTime_;

    bool isCurrentRunDead_; 
    unsigned long long currentRunDeathTime_;
    unsigned long long currentRunStartTime_; 
    
    unsigned long long totalEffectiveTimeNs_; 
    double totalEffectiveTimeSec_;            
    
    short currentProcessingRunID_; 
    bool suppressTimeoutWarnings_; 

    std::set<int> activeModuleIDs_;
    std::map<int, unsigned long long> moduleAliveTime_;
    unsigned long long prevEventTimeForEff_;

    std::map<int, long long> currentFileOffsets_;
    std::map<int, unsigned long long> baseTimeMap_;
    std::map<int, bool> hasSeenTimeHeaderMap_;

    std::vector<std::pair<unsigned long long, unsigned long long>> deadTimeRanges_;

    std::ofstream gainLogCsv_;
    std::ofstream rateLogCsv_;
    
    std::map<std::string, short> runPrefixToId_;
    std::vector<std::string> runIdToPrefix_;
    
    LUTItem ch_LUT[MAX_MODULES][MAX_SYS_CH];
    TH1F* hRawToT_LUT[MAX_MODULES][MAX_SYS_CH];
    TH1F* hGainCheck_LUT[MAX_MODULES][MAX_SYS_CH];

    std::map<std::pair<int, int>, ChConfig> detConfigMap_;
    std::map<std::pair<int, int>, TH1F*> hRawToTMap;
    std::map<std::pair<int, int>, TH1F*> gainCheckHists_;
    std::map<std::pair<int, int>, TH1F*> templateHists_;

    std::map<std::pair<int, int>, TGraph*> gainEvolutionGraphs_;
    std::map<std::pair<int, int>, TGraph*> rateEvolutionGraphs_;
    
    std::map<std::pair<int, int>, double> cumulativeShiftMap_; 
};

#endif