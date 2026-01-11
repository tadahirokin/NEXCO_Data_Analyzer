#ifndef DETECTOR_ANALYZER_H
#define DETECTOR_ANALYZER_H

#include <vector>
#include <string>
#include <map>
#include <set> // 追加
#include <deque>
#include <iostream>
#include <fstream>
#include <limits>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TSpectrum.h>

enum DetectorTypeID { ID_X1 = 0, ID_Y1 = 1, ID_X2 = 2, ID_Y2 = 3 };

struct Event {
    int type;       
    int strip;      
    int module;     
    unsigned long long eventTime_ns; 
    int tot;
    int uniqueKey;  
};

struct DetectorInfo {
    int detTypeID;
    int detInCh;
    std::string cachedName;
    int pIndex; 
};

class DetectorAnalyzer {
public:
    DetectorAnalyzer(int timeWindow_ns, const std::string& outputFileName);
    virtual ~DetectorAnalyzer();

    void loadConversionFactorsFromCSV(const std::string& csvFileName, const std::string& detName, int moduleID);
    void loadAndSortFiles(const std::string& directory, std::map<int, std::deque<std::string>>& fileQueues);
    void processBinaryFiles(std::map<int, std::deque<std::string>>& fileQueues);
    void setTimeLimitMinutes(double minutes);
    const std::map<int, std::vector<long long>>& getRawTotData() const;

private:
    void setupTree();
    
    std::pair<unsigned long long, bool> readEventsFromFile(
        const std::string& fileName, 
        std::vector<Event>& rawEvents, 
        long long& offset 
    );

    void processChunk(const std::vector<Event>& sortedEvents); 
    void analyzeGainShift(unsigned long long currentTime_ns);

    // 【追加】リアルタイムステータス表示用関数
    void printSearchStatus();

    std::string parseRunPrefix(const std::string& fileName);
    short getRunID(const std::string& prefix); 

    int timeWindow_ns_; 
    TFile* outputFile_;
    TTree* tree_;

    TH1F* hDeltaT_Nearest;
    TH1F* hDeltaT_n8;
    TH1F* hGlobalPIndexMap; 

    int     b_type;
    int     b_strip;
    long long b_time;
    int     b_tot;

    std::map<int, TH1F*> hRawToTMap;
    std::map<int, TH1F*> gainCheckHists_; 
    std::map<int, std::vector<long long>> rawTotDataVec_;
    std::map<int, DetectorInfo> uniqueKeyToDetInfo_;
    
    // 【追加】逆引きマップ: [検出器名][ストリップ番号] -> uniqueKey
    // 表を表示するために使います
    std::map<std::string, std::map<int, int>> nameStripToKey_;

    // 【追加】発見済みチャンネルの管理
    std::set<int> foundChannels_; 
    const int TOT_THRESHOLD_NS = 1000; // 1us = 1000ns以上のみ合格とする

    std::map<int, long long> currentFileOffsets_;

    unsigned long long current_base_time_ns_; 
    unsigned long long analysisDuration_ns_; 
    unsigned long long globalStartTimestamp_;

    const unsigned long long GAIN_ANALYSIS_WINDOW_NS = 600ULL * 1000000000ULL; 
    unsigned long long next_gain_analysis_time_ns_;
    std::map<int, double> refPeakPositions_;
    bool isFirstGainAnalysis_; 
    bool isCsvHeaderWritten_;
    std::ofstream gainLogCsv_; 
    std::ofstream rateLogCsv_; 

    std::vector<std::string> runIdToPrefix_;
    std::map<std::string, short> runPrefixToId_;
    
    std::map<int, long long> unknownKeyCounts_;
};

#endif