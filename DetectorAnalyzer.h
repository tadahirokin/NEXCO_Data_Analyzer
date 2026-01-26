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
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <chrono>
#include <iomanip>

// ROOT Headers
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
#include <TError.h>

// ID定義
const int ID_X1 = 0;
const int ID_Y1 = 1;
const int ID_X2 = 2;
const int ID_Y2 = 3;

// システム定数
const int MAX_MODULES = 16;   
const int MAX_SYS_CH = 512;
const int Y_CHANNEL_OFFSET = 64; 

// 閾値・設定
const unsigned long long GAIN_ANALYSIS_WINDOW_NS = 600000000000ULL; // 10分
const unsigned long long MODULE_DEATH_TIMEOUT_NS = 1000000000ULL;   // 1.0秒

// ヒストグラム設定 (20ns/bin)
const int MONITOR_HIST_BINS = 5000;      
const double MONITOR_HIST_MAX_TOT = 100000.0;
const double BIN_WIDTH_NS = 20.0; 

// ピーク探索・マッチング設定
const double PEAK_SEARCH_MIN_TOT = 10000.0; // 10us
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
        int globalIndex; 
    };

    DetectorAnalyzer(int timeWindow_ns, const std::string& outputFileName);
    ~DetectorAnalyzer();

    void setMuonAnalysisMode(bool enable);
    void loadConversionFactorsFromCSV(const std::string& csvFileName, const std::string& detName, int moduleID);
    void loadOfflineStripList(const std::string& csvFileName);
    void loadAndSortFiles(const std::string& directory, std::map<int, std::deque<std::string>>& fileQueues);
    void setTimeLimitMinutes(double minutes);
    void processBinaryFiles(std::map<int, std::deque<std::string>>& fileQueues);

    void printSearchStatus();
    bool isAnalysisStarted_; 

private:
    // --- コンストラクタ初期化リストで要求されている変数群 ---
    int timeWindow_ns_;
    TGraph* gT0Check;
    long t0GlobalCount_;
    unsigned long long firstT0Time_;
    unsigned long long lastTimeoutCheckTime_;
    unsigned long long currentRunDeathTime_;
    unsigned long long currentRunStartTime_;
    double totalEffectiveTimeSec_;
    unsigned long long prevEventTimeForEff_;

    // --- .cpp 1041行目付近で要求されているメイン配列 ---
    TH1F* hRawToT_Global[128];
    TH1F* hGainCheck_Global[128];

    // --- メソッド宣言 ---
    void printSurvivalReport(); // 278行目のエラーを解消
    void saveCurrentStatus(unsigned long long currentTime_ns);
    void setupTree();
    void generatePDFReport();
    void calculateEffectiveTime();
    void printLog(const std::string& msg);
    void printFileTail(const std::string& filePath, long long currentOffset);

    // --- 解析器の構成データ ---
    std::map<std::pair<int, int>, ChConfig> detConfigMap_; 
    std::map<std::string, std::map<int, std::pair<int, int>>> nameStripToMapKey_;

    // --- ヒストグラム & グラフ ---
    std::map<std::pair<int, int>, TH1F*> hRawToTMap;
    std::map<std::pair<int, int>, TH1F*> gainCheckHists_;
    std::map<std::pair<int, int>, TH1F*> templateHists_;
    std::map<std::pair<int, int>, TGraph*> gainEvolutionGraphs_;
    std::map<std::pair<int, int>, TGraph*> rateEvolutionGraphs_;

    // --- ゲイン追尾パラメータ ---
    std::map<std::pair<int, int>, double> initialPeakPosMap_;
    std::map<std::pair<int, int>, double> cumulativeShiftNsMap_;
    std::map<std::pair<int, int>, int> rangeMinMap_;
    std::map<std::pair<int, int>, int> rangeMaxMap_;

    // --- 状態監視 ---
    std::set<std::pair<int, int>> ignoredMonitorChannels_; 
    std::set<std::pair<int, int>> offlineStrips_;
    std::set<int> activeModuleIDs_;
    std::map<std::pair<int, int>, bool> aliveStatus_;
    std::map<std::pair<int, int>, std::set<int>> foundChannels_;
    std::map<int, unsigned long long> channelLastRefHitTime_;
    
    bool survivalReportShown_;
    bool dashboardShown_;
    bool isCurrentRunDead_;
    bool abortCurrentRun_;
    bool suppressTimeoutWarnings_;

    // --- 時間・統計管理 ---
    unsigned long long globalStartTimestamp_;
    unsigned long long currentEventTime_ns_;
    unsigned long long analysisDuration_ns_;
    unsigned long long totalEffectiveTimeNs_; 
    unsigned long long finalTotalTimeNs_;
    unsigned long long lastGlobalAnalysisTime_ns_;
    unsigned long long totalDataSize_;
    unsigned long long processedDataSize_;
    std::chrono::system_clock::time_point analysisStartTime_;
    
    bool isFirstGainAnalysis_;
    bool isTemplateCaptured_;
    bool isCsvHeaderWritten_;

    static const int COVELL_WINDOW_NS = 200;

    // --- Binary Sync & Decode ---
    unsigned long long findNextT0(const std::string& fileName, int modID, long long& offset);
    bool syncDataStream(std::ifstream& ifs, long long& skippedBytes);
    std::pair<unsigned long long, bool> readEventsFromFile(const std::string& fileName, int modID, std::vector<Event>& rawEvents, long long& offset);
    bool processChunk(const std::vector<Event>& sortedEvents);      
    unsigned long long getSafeTime(const std::map<int, unsigned long long>& lastTimes);
    std::string getRunSignature(const std::string& fileName);
    short getRunID(const std::string& prefix);
    std::string parseRunPrefix(const std::string& fullPath);

    // --- Gain Analysis ---
    bool analyzeGainShift(unsigned long long currentTime_ns);
    int findRightMostPeak(TH1F* hist); 
    void determineIntegrationRange(TH1F* hist, int peakBin, int& outMinBin, int& outMaxBin);
    int findBestShift(TH1F* hTarget, TH1F* hTemplate, int binMin, int binMax);
    double calculateResidual(TH1F* hTarget, TH1F* hTemplate, int shiftBins, int binMin, int binMax, double scale);
    double calculatePeakAreaCovell(TH1F* hist, double peakPos); 
    
    // --- ROOT I/O ---
    TFile* outputFile_;
    TTree* tree_;
    Char_t b_type;   
    Char_t b_strip;  
    Int_t  b_tot;    
    Long64_t b_time;   

    TH1F *hDeltaT_Neighbor;
    TH1F *hDeltaT_n8;
    TH1F *hGlobalPIndexMap;
    TH1F *hEventWidth_All;
    TH1F *hEventWidth_CondA;
    TH1F *hEventWidth_CondB;

    // --- Warp & Sync ---
    std::map<int, unsigned long long> lastT0_;     
    std::map<int, unsigned long long> lastRawT0_;  
    std::map<int, long long> moduleOffsets_;       
    std::map<int, bool> hasDataAligned_; 
    std::map<int, bool> hasFoundT0_;     
    std::map<int, long long> currentFileOffsets_;

    // --- Run ID ---
    std::map<std::string, short> runPrefixToId_;
    std::vector<std::string> runIdToPrefix_;
    short currentProcessingRunID_; 
    
    // --- LUT ---
    LUTItem ch_LUT[MAX_MODULES][MAX_SYS_CH];
    TH1F* hRawToT_LUT[MAX_MODULES][MAX_SYS_CH];
    TH1F* hGainCheck_LUT[MAX_MODULES][MAX_SYS_CH];

    // --- CSV Outputs ---
    std::ofstream gainLogCsv_;
    std::ofstream rateLogCsv_;

    // --- Muon Extraction ---
    void processEventExtraction(); 
    bool checkConditionA(const std::deque<Event>& buf, size_t startIdx, double& outWidth);
    bool checkConditionB(const std::deque<Event>& buf, size_t startIdx, double& outWidth);
    bool hasAdjacentPair(const std::vector<int>& strips);
    bool isAdjacentToOffline(int strip, int detTypeID);

    bool enableMuonAnalysis_;
    TFile* idealFile_;
    TTree* idealTree_;
    int    b_i_runID;
    int    b_i_eventID;
    int    b_i_nHits;
    double b_i_width;
    std::vector<int>    *b_i_type;
    std::vector<int>    *b_i_strip;
    std::vector<int>    *b_i_tot;
    std::vector<double> *b_i_time;
    std::deque<Event> analysisBuffer_;
    int globalEventID_Ideal_;

    std::map<int, unsigned long long> moduleAliveTime_;
    std::vector<std::pair<unsigned long long, unsigned long long>> deadTimeRanges_;
};

#endif