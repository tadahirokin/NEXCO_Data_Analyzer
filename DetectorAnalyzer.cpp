#include "DetectorAnalyzer.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cstdio>
#include <sys/stat.h>
#include <cstring>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include <regex>
#include <climits>
#include <TError.h>
#include <TROOT.h>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <thread>

#include <TLatex.h>
#include <TArrow.h>
#include <TBox.h>

namespace fs = std::filesystem;

// ----------------------------------------------------------------------------
// グローバル変数 & ログ関数
// ----------------------------------------------------------------------------
static DetectorAnalyzer *g_currentAnalyzer = nullptr;

// ----------------------------------------------------------------------------
// コンストラクタ
// ----------------------------------------------------------------------------
DetectorAnalyzer::DetectorAnalyzer(int timeWindow_ns, const std::string &outputFileName)
    : timeWindow_ns_(timeWindow_ns),
      outputFile_(new TFile(outputFileName.c_str(), "RECREATE")),
      tree_(nullptr), hDeltaT_Neighbor(nullptr), hDeltaT_n8(nullptr), hGlobalPIndexMap(nullptr),
      gT0Check(nullptr), t0GlobalCount_(0), firstT0Time_(0),
      analysisDuration_ns_(std::numeric_limits<unsigned long long>::max()),
      globalStartTimestamp_(0),
      currentEventTime_ns_(0),
      lastGlobalAnalysisTime_ns_(0),
      isFirstGainAnalysis_(true),
      isCsvHeaderWritten_(false),
      isTemplateCaptured_(false),
      totalDataSize_(0),
      processedDataSize_(0),
      dashboardShown_(false),
      lastTimeoutCheckTime_(0),
      isCurrentRunDead_(false),
      abortCurrentRun_(false),
      currentRunDeathTime_(0),
      currentRunStartTime_(0),
      totalEffectiveTimeNs_(0),
      totalEffectiveTimeSec_(0.0),
      currentProcessingRunID_(-1),
      suppressTimeoutWarnings_(false),
      prevEventTimeForEff_(0),
      enableMuonAnalysis_(false),
      idealFile_(nullptr), idealTree_(nullptr),
      hEventWidth_All(nullptr), hEventWidth_CondA(nullptr), hEventWidth_CondB(nullptr),
      b_i_type(nullptr), b_i_strip(nullptr), b_i_time(nullptr), b_i_corrected_tot(nullptr),
      globalEventID_Ideal_(0)
{
    g_currentAnalyzer = this;
    gROOT->SetBatch(kTRUE);

    analysisStartTime_ = std::chrono::system_clock::now();

    if (!outputFile_ || outputFile_->IsZombie())
    {
        throw std::runtime_error("[ERROR] Failed to create output ROOT file.");
    }

    for (int m = 0; m < MAX_MODULES; ++m)
    {
        for (int s = 0; s < MAX_SYS_CH; ++s)
        {
            ch_LUT[m][s].isValid = false;
            hRawToT_LUT[m][s] = nullptr;
            hGainCheck_LUT[m][s] = nullptr;
        }
    }

    for (int m : activeModuleIDs_)
    {
        lastT0_[m] = 0;
        hasDataAligned_[m] = false;
        hasFoundT0_[m] = false;

        // ★ これを追加
        lastRawT0_[m] = 0;
        moduleOffsets_[m] = 0;
    }

    for (int i = 0; i < 128; ++i)
    {
        currentScaleFactors_[i] = 1.0;   // 初期値は補正なし
        initialPeakPositions_[i] = 0.0;  // 未定
        hCorrectedToT_Global[i] = nullptr; 
    }

    setupTree();

    struct stat st = {};
    if (stat("./gainshift", &st) == -1)
    {
#ifdef _WIN32
        _mkdir("./gainshift");
#else
        mkdir("./gainshift", 0777);
#endif
    }
    gainLogCsv_.open("./gainshift/gain_monitoring.csv", std::ios::out);
    rateLogCsv_.open("./gainshift/rate_monitoring.csv", std::ios::out);

    printLog("Analyzer initialized. Output: " + outputFileName);
}

void DetectorAnalyzer::setMuonAnalysisMode(bool enable)
{
    enableMuonAnalysis_ = enable;
    if (enable)
    {
        printLog("[MODE] Muon Analysis ENABLED (Clustering Active, Slower)");
        
        idealFile_ = new TFile("MuonEvents_Ideal.root", "RECREATE");
        idealTree_ = new TTree("IdealMuons", "Extracted Muon Events");

        idealTree_->Branch("runID",   &b_i_runID,   "runID/I");
        idealTree_->Branch("eventID", &b_i_eventID, "eventID/I");
        idealTree_->Branch("nHits",   &b_i_nHits,   "nHits/I");
        idealTree_->Branch("width",   &b_i_width,   "width/D");

        b_i_type  = new std::vector<int>();
        b_i_strip = new std::vector<int>();
        // ★削除: b_i_tot = new std::vector<int>(); 
        b_i_time  = new std::vector<double>();
        b_i_corrected_tot = new std::vector<double>(); 

        idealTree_->Branch("type",  &b_i_type);
        idealTree_->Branch("strip", &b_i_strip);
        // ★削除: idealTree_->Branch("raw_tot", &b_i_tot);
        
        idealTree_->Branch("time",  &b_i_time);
        
        // 補正後ToTのみ保存
        idealTree_->Branch("corrected_tot", &b_i_corrected_tot);
    }
    else
    {
        printLog("[MODE] Fast Gain Check ONLY (Muon Extraction Disabled)");
    }
}

DetectorAnalyzer::~DetectorAnalyzer()
{
    printLog("Shutting down analyzer...");

    // 1. 最後のゲイン解析を実行
    if (currentEventTime_ns_ > lastGlobalAnalysisTime_ns_)
    {
        analyzeGainShift(currentEventTime_ns_);
    }

    // 2. 有効時間の最終計算とPDFレポートの生成
    calculateEffectiveTime();
    generatePDFReport();

    // 3. ROOTファイルの保存
    if (outputFile_ && outputFile_->IsOpen())
    {
        outputFile_->cd();
        if (tree_)
            tree_->Write();

        // 解析結果（ヒストグラム）を保存
        auto dirHist = outputFile_->mkdir("Histograms");
        if (dirHist)
        {
            dirHist->cd();
            for (auto const &[key, h] : hRawToTMap)
                if (h)
                    h->Write();
        }

        outputFile_->cd();
        auto dG = outputFile_->mkdir("GainHistory");
        if (dG)
        {
            dG->cd();
            for (auto const &[key, g] : gainEvolutionGraphs_)
                if (g && g->GetN() > 0)
                    g->Write();
        }

        outputFile_->Close();
        delete outputFile_;
        outputFile_ = nullptr;
    }

    // 4. ミュオン抽出結果の保存
    if (idealFile_ && idealFile_->IsOpen())
    {
        idealFile_->cd();
        if (idealTree_)
            idealTree_->Write();
        idealFile_->Close();
        delete idealFile_;
    }

    if (gainLogCsv_.is_open())
        gainLogCsv_.close();
    if (rateLogCsv_.is_open())
        rateLogCsv_.close();
}

static std::map<int, bool> g_allowNextT0Jump;

void DetectorAnalyzer::printSearchStatus()
{
    double liveTimeMin = (double)totalEffectiveTimeNs_ / 1.0e9 / 60.0;
    double progress = (analysisDuration_ns_ > 0 && analysisDuration_ns_ != ULLONG_MAX)
                          ? (double)totalEffectiveTimeNs_ / (double)analysisDuration_ns_ * 100.0
                          : 0.0;
    std::cout << "\r\033[K" << "Progress: [" << std::fixed << std::setprecision(1) << progress << "%] "
              << "Live Time: " << liveTimeMin << " min" << std::flush;
}

void DetectorAnalyzer::printSurvivalReport()
{
    std::cout << "\n"
              << std::string(70, '=') << "\n [SURVIVAL REPORT at 5.0 min]\n"
              << std::string(70, '-') << std::endl;
    std::string detNames[] = {"X1", "Y1", "X2", "Y2"};
    for (int t = 0; t < 4; ++t)
    {
        int alive = 0, total = 0;
        std::vector<int> silent;
        for (auto const &[key, h] : hRawToTMap)
        {
            if (ch_LUT[key.first][key.second].detTypeID == t)
            {
                total++;
                if (aliveStatus_[key])
                    alive++;
                else
                    silent.push_back(ch_LUT[key.first][key.second].strip);
            }
        }
        std::cout << " " << detNames[t] << ": " << alive << " / " << total << " Alive. ";
        if (!silent.empty())
        {
            std::cout << "Silent Strips: ";
            std::sort(silent.begin(), silent.end());
            for (size_t j = 0; j < silent.size(); ++j)
                std::cout << silent[j] << (j == silent.size() - 1 ? "" : ",");
        }
        std::cout << std::endl;
    }
    std::cout << std::string(70, '=') << "\n"
              << std::endl;
}


void DetectorAnalyzer::processEventExtraction()
{
    if (analysisBuffer_.empty()) return;

    while (analysisBuffer_.size() >= 8)
    {
        double w = 0;
        bool isA = checkConditionA(analysisBuffer_, 0, w);
        bool isB = isA ? false : checkConditionB(analysisBuffer_, 0, w);

        if (isA || isB)
        {
            b_i_runID = currentProcessingRunID_;
            b_i_eventID = globalEventID_Ideal_++;
            b_i_nHits = isA ? 8 : 7;
            b_i_width = w;

            b_i_type->clear();
            b_i_strip->clear();
            // ★削除: b_i_tot->clear();
            b_i_time->clear();
            b_i_corrected_tot->clear();

            for (int k = 0; k < b_i_nHits; ++k)
            {
                const auto &e = analysisBuffer_[k];
                b_i_type->push_back(ch_LUT[e.module][e.sysCh].detTypeID);
                b_i_strip->push_back(ch_LUT[e.module][e.sysCh].strip);
                
                // ★削除: b_i_tot->push_back(e.tot);
                
                b_i_time->push_back((double)e.eventTime_ns);
                
                // 補正後ToTを保存
                b_i_corrected_tot->push_back(e.correctedTot);
            }

            if (idealTree_) idealTree_->Fill();

            for (int k = 0; k < b_i_nHits; ++k) analysisBuffer_.pop_front();
        }
        else
        {
            analysisBuffer_.pop_front();
        }
    }
}

// 条件A: 全層で隣り合う2つを含む (8ヒット)
bool DetectorAnalyzer::checkConditionA(const std::deque<Event> &buf, size_t startIdx, double &outWidth)
{
    if (startIdx + 8 > buf.size())
        return false;
    outWidth = (double)(buf[startIdx + 7].eventTime_ns - buf[startIdx].eventTime_ns);

    std::map<int, std::vector<int>> sMap;
    for (size_t i = 0; i < 8; ++i)
    {
        const auto &e = buf[startIdx + i];
        int type = ch_LUT[e.module][e.sysCh].detTypeID;
        int strip = ch_LUT[e.module][e.sysCh].strip;
        sMap[type].push_back(strip);
    }
    // X1, Y1, X2, Y2 全てで隣接ペアを持つか
    for (int type = 0; type < 4; ++type)
    {
        if (!hasAdjacentPair(sMap[type]))
            return false;
    }
    return true;
}

// 条件B: Y1欠損考慮 (7ヒット)
bool DetectorAnalyzer::checkConditionB(const std::deque<Event> &buf, size_t startIdx, double &outWidth)
{
    if (startIdx + 7 > buf.size())
        return false;
    outWidth = (double)(buf[startIdx + 6].eventTime_ns - buf[startIdx].eventTime_ns);

    std::map<int, std::vector<int>> sMap;
    for (size_t i = 0; i < 7; ++i)
    {
        const auto &e = buf[startIdx + i];
        int type = ch_LUT[e.module][e.sysCh].detTypeID;
        int strip = ch_LUT[e.module][e.sysCh].strip;
        sMap[type].push_back(strip);
    }

    if (!hasAdjacentPair(sMap[ID_X1]))
        return false;
    if (!hasAdjacentPair(sMap[ID_X2]))
        return false;
    if (!hasAdjacentPair(sMap[ID_Y2]))
        return false;

    if (sMap[ID_Y1].size() != 1)
        return false;
    return isAdjacentToOffline(sMap[ID_Y1][0], ID_Y1);
}

bool DetectorAnalyzer::hasAdjacentPair(const std::vector<int> &strips)
{
    if (strips.size() < 2)
        return false;
    std::vector<int> s = strips;
    std::sort(s.begin(), s.end());
    for (size_t i = 0; i < s.size() - 1; ++i)
    {
        if (s[i + 1] - s[i] == 1)
            return true;
    }
    return false;
}

bool DetectorAnalyzer::isAdjacentToOffline(int strip, int detTypeID)
{
    if (offlineStrips_.count({detTypeID, strip - 1}))
        return true;
    if (offlineStrips_.count({detTypeID, strip + 1}))
        return true;
    return false;
}

// T0デコード (ビット演算一元化)
unsigned long long DetectorAnalyzer::decodeT0(const unsigned char *p)
{
    return ((unsigned long long)p[0] << 24 | (unsigned long long)p[1] << 16 |
            (unsigned long long)p[2] << 8 | (unsigned long long)p[3]) *
               1000000000ULL +
           ((unsigned long long)p[4] << 2 | (unsigned long long)(p[5] & 0xC0) >> 6) * 1000000ULL +
           ((unsigned long long)(p[5] & 0x3F) << 2 | (unsigned long long)(p[6] & 0xC0) >> 6);
}

// デッドタイム集計 (重複マージ)
unsigned long long DetectorAnalyzer::calculateTotalDeadTime(std::vector<std::pair<unsigned long long, unsigned long long>> &ranges)
{
    if (ranges.empty())
        return 0;
    std::sort(ranges.begin(), ranges.end());
    std::vector<std::pair<unsigned long long, unsigned long long>> merged;
    merged.push_back(ranges[0]);
    for (size_t i = 1; i < ranges.size(); ++i)
    {
        if (ranges[i].first < merged.back().second)
        {
            if (ranges[i].second > merged.back().second)
                merged.back().second = ranges[i].second;
        }
        else
            merged.push_back(ranges[i]);
    }
    unsigned long long total = 0;
    for (const auto &p : merged)
        total += (p.second - p.first);
    return total;
}

bool DetectorAnalyzer::performAlignment(const std::vector<unsigned char> &buf, size_t &idx, size_t maxLimit)
{
    for (; idx + 80 <= maxLimit; ++idx)
    {
        if (buf[idx] == 0x69 || buf[idx] == 0x6a)
        {
            bool match = true;
            for (int m = 1; m < 10; ++m)
            {
                if (buf[idx + m * 8] != 0x69 && buf[idx + m * 8] != 0x6a)
                {
                    match = false;
                    break;
                }
            }
            if (match)
                return true;
        }
    }
    return false;
}

int DetectorAnalyzer::findNextT0(const std::vector<unsigned char> &buf, size_t &idx, size_t maxLimit, unsigned long long &t0_out)
{
    for (; idx + 8 <= maxLimit; idx += 8)
    {
        if (buf[idx] == 0x69)
        {
            t0_out = decodeT0(&buf[idx + 1]);
            return 0; // FOUND
        }
        else if (buf[idx] != 0x6a)
            return -1; // LOST
    }
    return 1; // NOT_FOUND
}

long long DetectorAnalyzer::synchronizeStream(const std::string &fileName, int modID, long long startOffset, unsigned long long targetTime)
{
    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs.is_open())
        return startOffset;

    printLog(Form("[SYNC] Mod%d Start Sync from offset %lld", modID, startOffset));

    ifs.seekg(startOffset, std::ios::beg);
    const size_t SYNC_BUF_SIZE = 4 * 1024 * 1024;
    std::vector<unsigned char> buf(SYNC_BUF_SIZE);
    long long currentOffset = startOffset;
    bool aligned = false;

    while (ifs.good())
    {
        ifs.read(reinterpret_cast<char *>(buf.data()), SYNC_BUF_SIZE);
        size_t count = static_cast<size_t>(ifs.gcount());
        if (count < 80)
            break;
        size_t i = 0;
        if (!aligned)
        {
            if (performAlignment(buf, i, count))
                aligned = true;
            else
            {
                currentOffset += count;
                continue;
            }
        }
        unsigned long long t0_val = 0;
        bool bufferFinished = false;
        while (i + 8 <= count)
        {
            int res = findNextT0(buf, i, count, t0_val); // ここを findNextT0 に修正
            if (res == 0)
            {
                if (t0_val >= targetTime)
                    return currentOffset + (long long)i;
                i += 8;
            }
            else if (res == 1)
            {
                bufferFinished = true;
                break;
            }
            else
            {
                aligned = false;
                i++;
                break;
            }
        }
        if (bufferFinished)
            currentOffset += count;
        else
        {
            currentOffset += i;
            ifs.seekg(currentOffset, std::ios::beg);
        }
    }
    return currentOffset;
}


std::pair<unsigned long long, bool> DetectorAnalyzer::readEventsFromFile(const std::string &fileName, int modID, std::vector<Event> &rawEvents, long long &offset)
{
    // ========================================================================
    // プロセス1: アライメント (ファイル先頭時)
    // ========================================================================
    // 同一ラン内のファイル切り替えでも、ヘッダのビットズレ補正は必須。
    // targetTime=0 なので再同期（読み飛ばし）はせず、オフセット調整のみ行う。
    if (offset == 0)
    {
        offset = synchronizeStream(fileName, modID, 0, 0);
    }

    // ========================================================================
    // データ読み込み準備
    // ========================================================================
    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs.is_open()) return {lastT0_[modID], true};

    ifs.seekg(0, std::ios::end);
    long long fileSize = static_cast<long long>(ifs.tellg());
    
    if (offset >= fileSize) return {lastT0_[modID], true};

    const size_t CHUNK_SIZE = 64 * 1024 * 1024;
    ifs.seekg(offset, std::ios::beg);
    
    std::vector<unsigned char> buf(CHUNK_SIZE);
    ifs.read(reinterpret_cast<char *>(buf.data()), CHUNK_SIZE);
    size_t count = static_cast<size_t>(ifs.gcount());
    bool isEOF = (offset + static_cast<long long>(count) >= fileSize);

    if (count < 8)
    {
        offset += count;
        return {lastT0_[modID], isEOF};
    }

    size_t i = 0;
    unsigned long long currentBaseTime = lastT0_[modID];

    // ========================================================================
    // メインループ
    // ========================================================================
    while (i + 8 <= count)
    {
        unsigned char header = buf[i];

        // グリッドずれ検知
        if (header != 0x69 && header != 0x6a)
        {
            // ここでズレたら再アライメント
            offset = synchronizeStream(fileName, modID, offset + (long long)i + 1, 0);
            return {lastT0_[modID], isEOF};
        }

if (header == 0x69)
        {
            unsigned long long raw_t = decodeT0(&buf[i + 1]);
            
            // 初回（まだ一度もT0を読んでいない）の場合は比較対象がないのでスルー
            if (lastT0_[modID] != 0)
            {
                long long diff = (long long)raw_t - (long long)lastT0_[modID];

                // ★変更点: 1秒 (1e9 ns) 以上のタイムジャンプ検知
                // データ収集系のバグ（ファイル結合ミス）による「別ランの混入」とみなし、
                // このファイルの残りを全て破棄して強制終了する。
                if (std::abs(diff) > 1000000000LL)
                {
                    printLog("[WARN] Time Jump Detected in Mod " + std::to_string(modID) + 
                             " (Diff: " + std::to_string((double)diff/1e9) + " sec).");
                    printLog("       Last Valid: " + std::to_string(lastT0_[modID]) + 
                             " -> Found: " + std::to_string(raw_t));
                    printLog("       -> Truncating the rest of " + parseRunPrefix(fileName));

                    // ファイルポインタを末尾まで進めて「読み終わったこと」にする
                    offset = fileSize; 
                    
                    // ★重要: lastT0_ (現在時刻) は更新せずに戻る
                    return {lastT0_[modID], true}; 
                }
                
                // (通常のデッドタイム集計はそのまま残す)
                if (diff > 1100000LL)
                    deadTimeRanges_.push_back({lastT0_[modID], raw_t});
            }
            
            // 正常範囲なら時刻更新
            lastT0_[modID] = raw_t;
            currentBaseTime = raw_t;
            i += 8;
        }
        else if (header == 0x6a) // Data Packet (明示化)
        {
            unsigned char *p = &buf[i + 1];
            unsigned int tof = ((unsigned int)p[0] << 16 | (unsigned int)p[1] << 8 | (unsigned int)p[2]);
            unsigned int pw  = ((unsigned int)p[3] << 12 | (unsigned int)p[4] << 4 | (unsigned int)(p[5] & 0xF0) >> 4);
            
            // 警告対策の ", 0.0" を末尾に追加
            rawEvents.push_back({0, (int)p[6], modID, currentBaseTime + tof - pw, (int)pw, (int)p[6], 0.0});
            i += 8;
        }
    }

    offset += static_cast<long long>(i);
    return {lastT0_[modID], isEOF};
}

unsigned long long DetectorAnalyzer::getFileFirstT0(const std::string& filePath) {
    std::ifstream ifs(filePath, std::ios::binary);
    if (!ifs.is_open()) return 0;

    // スキャン範囲は広めに確保
    const size_t SCAN_SIZE = 8 * 1024 * 1024;
    std::vector<unsigned char> buf(SCAN_SIZE);
    ifs.read(reinterpret_cast<char*>(buf.data()), buf.size());
    size_t count = static_cast<size_t>(ifs.gcount());
    
    size_t idx = 0;
    while (idx + 160 <= count) { 
        // 1. 物理アライメント確認 (10パケット連続ヘッダ一致)
        if (performAlignment(buf, idx, count)) {
            std::vector<unsigned long long> foundT0s;
            size_t tempIdx = idx;

            // 2. その位置からT0(0x69)を10個抽出して検証
            while (foundT0s.size() < 10 && tempIdx + 8 <= count) {
                if (buf[tempIdx] == 0x69) {
                    foundT0s.push_back(decodeT0(&buf[tempIdx + 1]));
                } else if (buf[tempIdx] != 0x6a) {
                    break; // 未知のパケットがあればアライメントからやり直し
                }
                tempIdx += 8;
            }

            // 3. 10個のT0間隔が「2ms以内」に収まっているかチェック
            if (foundT0s.size() == 10) {
                bool isStable = true;
                for (size_t k = 0; k < foundT0s.size() - 1; ++k) {
                    long long diff = (long long)(foundT0s[k+1] - foundT0s[k]);
                    
                    // ★ 修正箇所：2ms (2,000,000ns) までの変動を許容
                    // 1ms周期に対し、1つパケットを落としても2ms弱になるはずなので、
                    // これにより微小な瞬きを許容しつつ、「年単位の狂い」は確実に弾きます。
                    // 0.1ms未満はノイズ（二重カウント等）の可能性があるため下限も設定。
                    if (diff < 100000LL || diff > 2000000LL) {
                        isStable = false;
                        break;
                    }
                }

                if (isStable) {
                    // 合格！このファイルを「物理的に連続した有効なデータ」とみなす
                    return foundT0s[0];
                }
            }
            idx++; // 条件を満たさなければ1バイトスライドして再探索
        } else {
            break; 
        }
    }
    return 0;
}

bool DetectorAnalyzer::processChunk(const std::vector<Event> &sortedEvents)
{
    if (sortedEvents.empty()) return true;

    if (prevEventTimeForEff_ == 0)
        prevEventTimeForEff_ = sortedEvents.front().eventTime_ns;

    for (const auto &e : sortedEvents)
    {
        // [時間管理] ラン開始時刻
        if (currentRunStartTime_ == 0 && e.tot >= 1000)
        {
            currentRunStartTime_ = e.eventTime_ns;
            if (lastGlobalAnalysisTime_ns_ == 0)
            {
                lastGlobalAnalysisTime_ns_ = e.eventTime_ns;
            }
            for (int mid : activeModuleIDs_)
                moduleAliveTime_[mid] = e.eventTime_ns;
        }

        currentEventTime_ns_ = e.eventTime_ns;
        moduleAliveTime_[e.module] = e.eventTime_ns;

        // ====================================================================
        // ★ここに戻しました: 10分経ったらその場で解析
        // ====================================================================
        if (lastGlobalAnalysisTime_ns_ > 0 && 
            e.eventTime_ns >= lastGlobalAnalysisTime_ns_ + GAIN_ANALYSIS_WINDOW_NS)
        {
            // 1. 解析実行
            // ここまでのデータを使って解析を行う
            analyzeGainShift(e.eventTime_ns);

            // 2. 基準時刻の更新
            // このイベントの時刻を「次の10分の開始地点」とする
            // これにより、もし直前に巨大なギャップがあっても、それはスキップされる
            lastGlobalAnalysisTime_ns_ = e.eventTime_ns;
        }

        // [ヒストグラム & Tree 充填]
        
// processChunk 内のループ修正

        // [ヒストグラム & Tree 充填]
        if (e.module >= 0 && e.module < MAX_MODULES && e.sysCh >= 0 && e.sysCh < MAX_SYS_CH)
        {
            // まずGlobal Indexを計算
            int typeID = ch_LUT[e.module][e.sysCh].detTypeID;
            int strip = ch_LUT[e.module][e.sysCh].strip;
            int gIdx = (typeID * 32) + (strip - 1);

            TH1F *h = hRawToT_LUT[e.module][e.sysCh];
            if (h != nullptr)
            {
                h->Fill((double)e.tot);

                // ゲイン監視用
                if (hGainCheck_LUT[e.module][e.sysCh])
                {
                    hGainCheck_LUT[e.module][e.sysCh]->Fill((double)e.tot);
                }

                aliveStatus_[{e.module, e.sysCh}] = true;

                // ★追加: ゲイン補正計算とヒストグラム充填
                double correctedVal = (double)e.tot;
                if (gIdx >= 0 && gIdx < 128)
                {
                    double scale = currentScaleFactors_[gIdx];
                    correctedVal = (double)e.tot * scale;
                    b_corrected_tot = correctedVal; // Tree用

                    if (hCorrectedToT_Global[gIdx])
                    {
                        fillFluxConserving(hCorrectedToT_Global[gIdx], (double)e.tot, scale);
                    }
                }

                if (hGlobalPIndexMap)
                {
                    hGlobalPIndexMap->Fill(gIdx);
                }

                if (tree_)
                {
                    b_type = static_cast<Char_t>(typeID);
                    b_strip = static_cast<Char_t>(strip);
                    b_tot = static_cast<Int_t>(e.tot);
                    b_time = static_cast<Long64_t>(e.eventTime_ns);
                    tree_->Fill();
                }
                
                // ★Muon解析用に補正値をセットするためにイベントをコピー更新する必要があるが
                // ここではループ変数のeがconst参照なので、バッファに入れるときにセットする
            }
        }

        if (enableMuonAnalysis_)
        {
            Event e_copy = e;
            // ★追加: 再計算した補正値をセット
            int typeID = ch_LUT[e.module][e.sysCh].detTypeID;
            int strip = ch_LUT[e.module][e.sysCh].strip;
            int gIdx = (typeID * 32) + (strip - 1);
            if (gIdx >= 0 && gIdx < 128) {
                e_copy.correctedTot = (double)e.tot * currentScaleFactors_[gIdx];
            } else {
                e_copy.correctedTot = (double)e.tot;
            }
            analysisBuffer_.push_back(e_copy);
        }
        prevEventTimeForEff_ = e.eventTime_ns;
    }

    // [Delta T]
    size_t N = sortedEvents.size();
    if (hDeltaT_Neighbor && N >= 2)
    {
        for (size_t i = 0; i < N - 1; ++i)
        {
            long long diff = (long long)(sortedEvents[i + 1].eventTime_ns - sortedEvents[i].eventTime_ns);
            if (diff >= 0) hDeltaT_Neighbor->Fill(diff);
        }
    }
    if (hDeltaT_n8 && N >= 8)
    {
        for (size_t i = 0; i <= N - 8; ++i)
        {
            long long diff = (long long)(sortedEvents[i + 7].eventTime_ns - sortedEvents[i].eventTime_ns);
            if (diff >= 0) hDeltaT_n8->Fill(diff);
        }
    }

    if (enableMuonAnalysis_)
    {
        processEventExtraction();
    }

    return true;
}

void DetectorAnalyzer::processBinaryFiles(std::map<int, std::deque<std::string>> &fileQueues)
{
    const size_t PER_MOD_CAP = 10000000;
    std::map<int, std::deque<Event>> moduleBuffers;
    std::map<int, unsigned long long> moduleLastTime;

    auto lastUIDraw = std::chrono::system_clock::now();

    unsigned long long accumulatedLiveTime_ns = 0;
    unsigned long long currentRunStartSafeTime = 0;
    unsigned long long lastLoopSafeTime = 0;
    unsigned long long currentRunGapTime_ns = 0;
    bool isFirstSafeTime = true;

    activeModuleIDs_.clear();
    for (auto const &[m, q] : fileQueues)
        if (!q.empty()) activeModuleIDs_.insert(m);

    for (int m : activeModuleIDs_)
    {
        lastT0_[m] = 0;
        hasDataAligned_[m] = false;
        hasFoundT0_[m] = false;
        moduleLastTime[m] = 0;
        currentFileOffsets_[m] = 0;
        lastRawT0_[m] = 0;
        moduleOffsets_[m] = 0;
    }

    std::string currentRunSignature = "";
    for (auto const &[m, q] : fileQueues)
        if (!q.empty())
        {
            currentRunSignature = getRunSignature(q.front());
            break;
        }

    // ========================================================================
    // ★初期同期 (Program Start)
    // ========================================================================
    printLog("[SYNC] Starting Initial Synchronization Sequence...");
    unsigned long long maxInitialT0 = 0;
    std::map<int, unsigned long long> initialT0s;

    for (int m : activeModuleIDs_)
    {
        if (fileQueues[m].empty()) continue;
        
        std::string fPath = fileQueues[m].front();
        
        // 1. アライメント (8バイト境界の確立)
        long long alignedOffset = synchronizeStream(fPath, m, 0, 0);
        currentFileOffsets_[m] = alignedOffset;

        // 2. 信頼できる最初のT0を探索 (Validation付き)
        std::ifstream ifs(fPath, std::ios::binary);
        if (ifs.is_open())
        {
            ifs.seekg(alignedOffset, std::ios::beg);
            std::vector<unsigned char> buf(2 * 1024 * 1024); // 2MB探索
            ifs.read(reinterpret_cast<char*>(buf.data()), buf.size());
            size_t count = ifs.gcount();
            size_t idx = 0;
            
            unsigned long long foundTime = 0;
            bool validT0Found = false;

            while (idx + 16 <= count) // 次のパケットも見たいので余裕を持つ
            {
                if (buf[idx] == 0x69)
                {
                    unsigned long long t0 = decodeT0(&buf[idx + 1]);
                    
                    // ★妥当性チェック: このT0がノイズでないか、少し先を見て確認
                    // 次のT0が「あまりにも遠い未来」や「過去」でないか簡易チェック
                    // (厳密すぎると読めなくなるので、あくまでデータ化け対策)
                    if (t0 > 0 && t0 < 18446744073709551615ULL) 
                    {
                        foundTime = t0;
                        validT0Found = true;
                        currentFileOffsets_[m] = alignedOffset + idx; // ここを開始位置にする
                        break;
                    }
                }
                // グリッドを進める（0x6aならスキップ、それ以外なら1バイト進んで再アライメント的な動き）
                if (buf[idx] == 0x69 || buf[idx] == 0x6a) idx += 8;
                else idx++;
            }

            if (validT0Found)
            {
                initialT0s[m] = foundTime;
                if (foundTime > maxInitialT0) maxInitialT0 = foundTime;
            }
            else
            {
                printLog("[WARN] No valid T0 found in header of M" + std::to_string(m));
                initialT0s[m] = 0;
            }
        }
    }

    // 3. T0同期 (遅れている方を読み飛ばす)
    printLog("[SYNC] Target Start Time: " + std::to_string(maxInitialT0));
    for (int m : activeModuleIDs_)
    {
        if (fileQueues[m].empty()) continue;

        if (initialT0s[m] > 0 && initialT0s[m] < maxInitialT0)
        {
            long long syncedOffset = synchronizeStream(fileQueues[m].front(), m, currentFileOffsets_[m], maxInitialT0);
            currentFileOffsets_[m] = syncedOffset;
        }
        
        // ★初期化時はlastT0にセット
        lastT0_[m] = maxInitialT0;
        hasFoundT0_[m] = true;
        hasDataAligned_[m] = true;
        moduleLastTime[m] = maxInitialT0;
    }
    
    lastGlobalAnalysisTime_ns_ = maxInitialT0;
    currentRunStartSafeTime = maxInitialT0;
    printLog("[SYNC] Synchronization Complete.");


    while (true)
    {
        // ... (ラン切り替え判定ロジック) ...
        bool allBuffersEmpty = true;
        bool allCurrentFilesDone = true;

        for (int m : activeModuleIDs_)
        {
            if (!moduleBuffers[m].empty()) allBuffersEmpty = false;
            if (!fileQueues[m].empty())
            {
                if (getRunSignature(fileQueues[m].front()) == currentRunSignature)
                    allCurrentFilesDone = false;
            }
        }

        if (allCurrentFilesDone && !allBuffersEmpty)
        {
            printLog("[DEBUG] Force flushing remaining buffers...");
            for (int m : activeModuleIDs_) if (!moduleBuffers[m].empty()) moduleBuffers[m].clear();
            allBuffersEmpty = true;
        }

        if (allBuffersEmpty)
        {
            std::string nextRunSignature = "";
            bool hasMoreFiles = false;
            for (int m : activeModuleIDs_)
            {
                if (!fileQueues[m].empty())
                {
                    hasMoreFiles = true;
                    std::string sig = getRunSignature(fileQueues[m].front());
                    if (nextRunSignature.empty()) nextRunSignature = sig;
                }
            }

            if (!hasMoreFiles) break;

            if (nextRunSignature != currentRunSignature)
            {
                printLog("[DEBUG] Run Switch Triggered! -> " + nextRunSignature);

                if (lastLoopSafeTime > lastGlobalAnalysisTime_ns_)
                {
                    unsigned long long tailDur = lastLoopSafeTime - lastGlobalAnalysisTime_ns_;
                    if (tailDur >= GAIN_ANALYSIS_WINDOW_NS * 0.9) analyzeGainShift(lastLoopSafeTime);
                }
                for (auto const &[key, h] : gainCheckHists_) if (h) h->Reset();
                printLog("[DEBUG] Gain Histograms flushed.");

                if (currentRunStartSafeTime > 0 && lastLoopSafeTime >= currentRunStartSafeTime)
                {
                    unsigned long long liveTime = (lastLoopSafeTime > currentRunStartSafeTime) ? (lastLoopSafeTime - currentRunStartSafeTime - currentRunGapTime_ns) : 0;
                    accumulatedLiveTime_ns += liveTime;
                    printLog(" [RUN RESULT] Live Time: " + std::to_string((double)liveTime / 1e9) + " sec");
                }

                currentRunSignature = nextRunSignature;

                // ------------------------------------------------------------
                // ★ラン切り替え時の同期 (Run-Switch Sync Ritual)
                // ------------------------------------------------------------
                printLog("[SYNC] Resyncing for New Run...");
                
                for (int m : activeModuleIDs_) {
                    lastT0_[m] = 0; 
                    hasFoundT0_[m] = false;
                    hasDataAligned_[m] = false;
                    currentFileOffsets_[m] = 0;
                    moduleLastTime[m] = 0;
                }

                maxInitialT0 = 0;
                initialT0s.clear();

                // Step 1 & 2: アライメント + 信頼できるT0探索
                for (int m : activeModuleIDs_)
                {
                    if (fileQueues[m].empty()) continue;
                    
                    long long alignedOffset = synchronizeStream(fileQueues[m].front(), m, 0, 0);
                    currentFileOffsets_[m] = alignedOffset;

                    std::ifstream ifs(fileQueues[m].front(), std::ios::binary);
                    if (ifs.is_open())
                    {
                        ifs.seekg(alignedOffset, std::ios::beg);
                        std::vector<unsigned char> buf(2 * 1024 * 1024);
                        ifs.read(reinterpret_cast<char*>(buf.data()), buf.size());
                        size_t count = ifs.gcount();
                        size_t idx = 0;
                        
                        unsigned long long foundTime = 0;
                        bool validT0Found = false;

                        while (idx + 16 <= count)
                        {
                            if (buf[idx] == 0x69)
                            {
                                unsigned long long t0 = decodeT0(&buf[idx + 1]);
                                // 簡易妥当性チェック
                                if (t0 > 0 && t0 < 18446744073709551615ULL)
                                {
                                    foundTime = t0;
                                    validT0Found = true;
                                    currentFileOffsets_[m] = alignedOffset + idx;
                                    break;
                                }
                            }
                            if (buf[idx] == 0x69 || buf[idx] == 0x6a) idx += 8;
                            else idx++;
                        }

                        if (validT0Found)
                        {
                            initialT0s[m] = foundTime;
                            if (foundTime > maxInitialT0) maxInitialT0 = foundTime;
                        }
                        else
                        {
                            initialT0s[m] = 0;
                        }
                    }
                }

                // Step 3: 同期
                printLog("[SYNC] New Run Target: " + std::to_string(maxInitialT0));
                for (int m : activeModuleIDs_)
                {
                    if (fileQueues[m].empty()) continue;

                    if (initialT0s[m] > 0 && initialT0s[m] < maxInitialT0)
                    {
                        long long syncedOffset = synchronizeStream(fileQueues[m].front(), m, currentFileOffsets_[m], maxInitialT0);
                        currentFileOffsets_[m] = syncedOffset;
                    }
                    
                    // ★重要: ラン切り替え後なので、同期した絶対時刻をセットする
                    lastT0_[m] = maxInitialT0;
                    hasFoundT0_[m] = true;
                    hasDataAligned_[m] = true;
                    moduleLastTime[m] = maxInitialT0;
                }

                currentRunStartSafeTime = maxInitialT0;
                lastGlobalAnalysisTime_ns_ = 0; // ギャップ除外
                lastLoopSafeTime = maxInitialT0;
                currentRunGapTime_ns = 0;
                isFirstSafeTime = false;

                continue;
            }
        }

        // ====================================================================
        // データ読み込み
        // ====================================================================
        int lagMod = -1;
        unsigned long long minT = 0xFFFFFFFFFFFFFFFFULL;
        
        for (int m : activeModuleIDs_)
        {
            if (fileQueues[m].empty() && moduleBuffers[m].empty()) continue;
            if (moduleLastTime[m] < minT)
            {
                minT = moduleLastTime[m];
                lagMod = m;
            }
        }

        if (lagMod == -1) break;

        if (lagMod != -1 && moduleBuffers[lagMod].size() < PER_MOD_CAP && !fileQueues[lagMod].empty())
        {
            std::string fileSig = getRunSignature(fileQueues[lagMod].front());
            if (fileSig == currentRunSignature)
            {
                std::vector<Event> temp;
                auto res = readEventsFromFile(fileQueues[lagMod].front(), lagMod, temp, currentFileOffsets_[lagMod]);
                if (!temp.empty())
                {
                    moduleBuffers[lagMod].insert(moduleBuffers[lagMod].end(), std::make_move_iterator(temp.begin()), std::make_move_iterator(temp.end()));
                    processedDataSize_ += (64 * 1024 * 1024);
                }
                moduleLastTime[lagMod] = res.first;
                
                if (res.second)
                {
                    printLog("[FILE DONE] " + fileQueues[lagMod].front());
                    fileQueues[lagMod].pop_front();
                    currentFileOffsets_[lagMod] = 0;
                }
            }
            else
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                continue; 
            }
        }

        // --- 同期時刻 (SafeTime) ---
        unsigned long long safeTime = getSafeTime(moduleLastTime);


// ---------------------------------------------------------
        // processBinaryFiles 内の if (safeTime > 0) ブロック
        // ---------------------------------------------------------
        if (safeTime > 0)
        {
            // [1. 時間管理変数の更新]
            if (isFirstSafeTime || safeTime < lastLoopSafeTime)
            {
                currentRunStartSafeTime = safeTime;
                isFirstSafeTime = false;
                lastGlobalAnalysisTime_ns_ = currentRunStartSafeTime;
            }
            if (lastGlobalAnalysisTime_ns_ == 0) lastGlobalAnalysisTime_ns_ = safeTime;

            // [2. 経過時間の計算と制限チェック]
            unsigned long long currentRunWallTime = (safeTime >= currentRunStartSafeTime) ? (safeTime - currentRunStartSafeTime) : 0;
            unsigned long long currentRunLiveTime = (currentRunWallTime > currentRunGapTime_ns) ? (currentRunWallTime - currentRunGapTime_ns) : 0;

            finalTotalTimeNs_ = accumulatedLiveTime_ns + currentRunLiveTime;
            lastLoopSafeTime = safeTime;

            if (analysisDuration_ns_ != std::numeric_limits<unsigned long long>::max())
            {
                if (finalTotalTimeNs_ >= analysisDuration_ns_) break;
            }

            // [3. データをバッファから取り出して処理]
            std::vector<Event> chunk;
            for (int m : activeModuleIDs_)
            {
                while (!moduleBuffers[m].empty() && moduleBuffers[m].front().eventTime_ns <= safeTime)
                {
                    chunk.push_back(std::move(moduleBuffers[m].front()));
                    moduleBuffers[m].pop_front();
                }
            }
            
            if (!chunk.empty())
            {
                std::sort(chunk.begin(), chunk.end(), [](const Event &a, const Event &b)
                          { return a.eventTime_ns < b.eventTime_ns; });
                
                // 作業員にデータを渡す
                // ★解析の実行判定は、この processChunk の中で行われるようになったので、
                //   ここで analyzeGainShift を呼ぶ必要はもうありません。
                processChunk(chunk); 
            }

            // ★以前ここにあった「analyzeGainShift呼び出しブロック」は完全に削除！
        }

        // --- UI更新 (前回と同じなので省略) ---
        static std::map<int, std::string> lastTrackedFile;
        auto now = std::chrono::system_clock::now();
        bool fileSwitched = false;
        std::string switchReason = "";

        for (int m : activeModuleIDs_)
        {
            std::string currentName = "";
            if (!fileQueues[m].empty()) currentName = std::filesystem::path(fileQueues[m].front()).filename().string();
            else if (!moduleBuffers[m].empty()) currentName = "(Buffering Last)";
            else currentName = "(Finished)";

            if (lastTrackedFile[m].empty()) lastTrackedFile[m] = currentName;
            if (lastTrackedFile[m] != currentName) {
                fileSwitched = true;
                switchReason = "Mod " + std::to_string(m) + " Changed";
            }
        }

        if (fileSwitched)
        {
            std::cout << "\n\033[K" << "---------------------------------------------------------------------------" << std::endl;
            std::cout << "\033[K" << "[FILE SWITCH] " << switchReason << " (SysTime: "
                      << std::fixed << std::setprecision(1) << (double)safeTime / 1e9 / 60.0 << " m)" << std::endl;

            for (int m : activeModuleIDs_)
            {
                std::string fPath = (!fileQueues[m].empty()) ? fileQueues[m].front() : "Done";
                std::string fName = std::filesystem::path(fPath).filename().string();
                if (fPath == "Done") fName = "(Finished)";

                long long currentAddr = currentFileOffsets_[m];
                long long fSize = (fPath != "Done" && std::filesystem::exists(fPath)) ? std::filesystem::file_size(fPath) : 1;
                double pct = (fPath == "Done") ? 100.0 : (double)currentAddr / (double)fSize * 100.0;
                double bufRatio = (double)moduleBuffers[m].size() / PER_MOD_CAP;
                
                int filledB = (int)(std::min(1.0, bufRatio) * 10);
                std::string barB = "["; for (int i = 0; i < 10; ++i) barB += (i < filledB ? '#' : '.'); barB += "]";
                int filledF = (int)(std::min(100.0, pct) / 100.0 * 10);
                std::string barF = "["; for (int i = 0; i < 10; ++i) barF += (i < filledF ? '=' : '.'); barF += "]";

                if (fName.length() > 20) fName = "..." + fName.substr(fName.length() - 17);
                double tMin = (double)lastT0_[m] / 1e9 / 60.0;

                std::cout << "\033[K"
                          << " M" << m << " | " << std::left << std::setw(20) << fName
                          << " | F:" << barF << std::right << std::setw(5) << std::fixed << std::setprecision(1) << pct << "%"
                          << " | T:" << std::setw(9) << std::fixed << std::setprecision(1) << tMin << "m"
                          << " | B:" << barB << std::endl;

                if (!fileQueues[m].empty()) lastTrackedFile[m] = std::filesystem::path(fileQueues[m].front()).filename().string();
                else lastTrackedFile[m] = "(Finished)";
            }
            std::cout << "\033[K" << "---------------------------------------------------------------------------" << std::endl;
            lastUIDraw = std::chrono::time_point<std::chrono::system_clock>::min();
        }

        if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastUIDraw).count() > 200)
        {
            lastUIDraw = now;
            std::stringstream ss;
            int lineCount = 0;
            for (int m : activeModuleIDs_)
            {
                if (lineCount > 0) ss << "\n";
                std::string fPath = (!fileQueues[m].empty()) ? fileQueues[m].front() : "Done";
                std::string fName = std::filesystem::path(fPath).filename().string();
                if (fPath == "Done") fName = "(Finished)";
                long long currentAddr = currentFileOffsets_[m];
                long long fSize = (fPath != "Done" && std::filesystem::exists(fPath)) ? std::filesystem::file_size(fPath) : 1;
                double pct = (fPath == "Done") ? 100.0 : (double)currentAddr / (double)fSize * 100.0;
                double bufRatio = (double)moduleBuffers[m].size() / PER_MOD_CAP;
                int filledF = (int)(std::min(100.0, pct) / 100.0 * 10);
                std::string barF = "["; for (int i = 0; i < 10; ++i) barF += (i < filledF ? '=' : '.'); barF += "]";
                int filledB = (int)(std::min(1.0, bufRatio) * 10);
                std::string barB = "["; for (int i = 0; i < 10; ++i) barB += (i < filledB ? '#' : '.'); barB += "]";
                if (fName.length() > 20) fName = "..." + fName.substr(fName.length() - 17);
                double tMin = (double)lastT0_[m] / 1e9 / 60.0;
                ss << "\033[K" << " M" << m << " | " << std::left << std::setw(20) << fName 
                   << " | F:" << barF << std::right << std::setw(5) << std::fixed << std::setprecision(1) << pct << "%" 
                   << " | T:" << std::setw(9) << std::fixed << std::setprecision(1) << tMin << "m" 
                   << " | B:" << barB;
                lineCount++;
            }
            std::cout << ss.str() << std::flush;
            if (lineCount > 1) std::cout << "\033[" << (lineCount - 1) << "A";
            std::cout << "\r" << std::flush;
        }
    }
    std::cout << "\033[" << (activeModuleIDs_.size() + 1) << "B" << std::endl;
}

// ★修正: 戻り値を bool に変更 (ヘッダと一致させる)
bool DetectorAnalyzer::analyzeGainShift(unsigned long long endTime_ns)
{
    // 1. 経過時間の計算
    double durationSec = 0.0;
    if (lastGlobalAnalysisTime_ns_ > 0 && endTime_ns > lastGlobalAnalysisTime_ns_)
    {
        durationSec = (double)(endTime_ns - lastGlobalAnalysisTime_ns_) / 1.0e9;
    }

    // 2. 統計不足判定 (90%未満なら捨てる)
    double targetWindowSec = (double)GAIN_ANALYSIS_WINDOW_NS / 1.0e9;
    if (durationSec < (targetWindowSec * 0.9))
    {
        for (auto const &[key, h] : gainCheckHists_)
        {
            if (h) h->Reset(); 
        }
        return false; // ★修正: falseを返す
    }

    double plotTimeMin = (double)lastGlobalAnalysisTime_ns_ / 60.0e9;

    // CSVヘッダ出力
    if (!isCsvHeaderWritten_)
    {
        std::string detNames[] = {"X1", "Y1", "X2", "Y2"};
        gainLogCsv_ << "Time_ns";
        rateLogCsv_ << "Time_ns";
        for (int type = 0; type < 4; ++type)
        {
            for (int s = 1; s <= 32; ++s)
            {
                std::string header = "," + detNames[type] + "_" + std::to_string(s);
                gainLogCsv_ << header;
                rateLogCsv_ << header;
            }
        }
        gainLogCsv_ << "\n";
        rateLogCsv_ << "\n";
        isCsvHeaderWritten_ = true;
    }

    std::vector<double> currentStepPeaks(128, -1.0);
    std::vector<double> currentStepRates(128, 0.0);

    for (auto const &[histKey, hist] : gainCheckHists_)
    {
        if (!hist) continue;
        if (hist->GetEntries() < 10)
        {
            hist->Reset();
            continue;
        }

        // 新関数呼び出し
        double peakNs = 0.0;
        double totalCounts = 0.0;
        bool found = findPeakAndIntegrate(hist, peakNs, totalCounts);

        if (found && peakNs > 0.0)
        {
            int idx = detConfigMap_[histKey].detTypeID * 32 + (detConfigMap_[histKey].strip - 1);
            
            if (idx >= 0 && idx < 128) currentStepPeaks[idx] = peakNs;

            // グラフ更新
            if (!gainEvolutionGraphs_[histKey])
            {
                gainEvolutionGraphs_[histKey] = new TGraph();
                gainEvolutionGraphs_[histKey]->SetTitle(Form("Peak Trend: %s_%d", detConfigMap_[histKey].cachedName.c_str(), detConfigMap_[histKey].strip));
                gainEvolutionGraphs_[histKey]->GetXaxis()->SetTitle("Time [min]");
                gainEvolutionGraphs_[histKey]->GetYaxis()->SetTitle("Peak Position [ns]");
            }
            gainEvolutionGraphs_[histKey]->SetPoint(gainEvolutionGraphs_[histKey]->GetN(), plotTimeMin, peakNs);

            if (!gainOverlayGraphs_[histKey])
            {
                gainOverlayGraphs_[histKey] = new TGraph();
                gainOverlayGraphs_[histKey]->SetTitle(Form("Peak Drift: %s_%d", detConfigMap_[histKey].cachedName.c_str(), detConfigMap_[histKey].strip));
                gainOverlayGraphs_[histKey]->GetXaxis()->SetTitle("ToT [ns]");
                gainOverlayGraphs_[histKey]->GetYaxis()->SetTitle("Time [min]");
            }
            gainOverlayGraphs_[histKey]->SetPoint(gainOverlayGraphs_[histKey]->GetN(), peakNs, plotTimeMin);

            if (idx >= 0 && idx < 128)
            {
                currentStepRates[idx] = (totalCounts / durationSec) * 60.0;
            }

            if (!rateEvolutionGraphs_[histKey])
            {
                rateEvolutionGraphs_[histKey] = new TGraph();
                rateEvolutionGraphs_[histKey]->SetTitle(Form("Rate Trend: %s_%d", detConfigMap_[histKey].cachedName.c_str(), detConfigMap_[histKey].strip));
                rateEvolutionGraphs_[histKey]->GetXaxis()->SetTitle("Time [min]");
                rateEvolutionGraphs_[histKey]->GetYaxis()->SetTitle("Rate [CPM]");
            }
            rateEvolutionGraphs_[histKey]->SetPoint(rateEvolutionGraphs_[histKey]->GetN(), plotTimeMin, currentStepRates[idx]);
        }
        
        hist->Reset();
    }

    // CSV出力
    gainLogCsv_ << lastGlobalAnalysisTime_ns_;
    rateLogCsv_ << lastGlobalAnalysisTime_ns_;

    for (int i = 0; i < 128; ++i)
    {
        if (currentStepPeaks[i] > 0.0)
            gainLogCsv_ << "," << std::fixed << std::setprecision(2) << currentStepPeaks[i];
        else
            gainLogCsv_ << ",ND";
        rateLogCsv_ << "," << std::fixed << std::setprecision(2) << currentStepRates[i];
    }
    gainLogCsv_ << "\n";
    rateLogCsv_ << "\n";
    gainLogCsv_.flush();
    rateLogCsv_.flush();

    return true; // ★修正: trueを返す
}

// ----------------------------------------------------------------------------
// 結果表示 (LiveTime表示の更新)
// ----------------------------------------------------------------------------
void DetectorAnalyzer::calculateEffectiveTime()
{
    double liveTimeSec = (double)finalTotalTimeNs_ / 1.0e9;

    // 設定時間
    double targetSec = (analysisDuration_ns_ != std::numeric_limits<unsigned long long>::max())
                           ? (double)analysisDuration_ns_ / 1.0e9
                           : 0.0;

    std::cout << "\n========================================================" << std::endl;
    std::cout << " [RESULT] Total System-wide Live Time: " << std::fixed << std::setprecision(2) << liveTimeSec << " sec" << std::endl;
    std::cout << "          (Accumulated from synchronized T0 intervals)" << std::endl;

    if (targetSec > 0)
    {
        double ratio = (liveTimeSec / targetSec) * 100.0;
        std::cout << "          Target Completion Ratio: " << std::fixed << std::setprecision(2) << ratio << " % (Limit: " << targetSec << "s)" << std::endl;
    }
    std::cout << "========================================================\n"
              << std::endl;
}


bool DetectorAnalyzer::findPeakAndIntegrate(TH1F* h, double& outCentroid, double& outCount)
{
    if (!h) return false;

    int nBins = h->GetNbinsX();
    int silenceThreshold = (int)(2000.0 / BIN_WIDTH_NS); // 2µs
    
    // ★修正: 信号本体(Core)の閾値を 0.2µs (200ns) に変更
    // ノイズは0.2µs以下なので、それより大きければ信号とみなす
    double minCoreSignalWidthNs = 200.0; 

    // ペデスタル領域(10µs以下)には入らない
    int minSearchBin = (int)(PEAK_SEARCH_MIN_TOT / BIN_WIDTH_NS); 

    int searchHead = nBins;

    while (searchHead > minSearchBin)
    {
        // 1. 信号の「右端」を探す
        int signalRightBin = -1;
        for (int i = searchHead; i >= minSearchBin; --i)
        {
            if (h->GetBinContent(i) > 0)
            {
                signalRightBin = i;
                break;
            }
        }
        if (signalRightBin == -1) return false;

        // 2. 左側へスキャンし、2µsの空白を探す
        int signalLeftBin = 1; 
        int integLeftBin = 1;  
        int zeroCount = 0;
        bool gapFound = false;

        for (int i = signalRightBin; i >= minSearchBin; --i)
        {
            if (h->GetBinContent(i) <= 0) zeroCount++;
            else zeroCount = 0;

            if (zeroCount >= silenceThreshold)
            {
                integLeftBin = i; 
                // 信号本体の左端 (空白が終わった直後)
                signalLeftBin = i + zeroCount; 
                gapFound = true;
                break;
            }
        }
        
        // 探索下限(10us)まで空白が見つからなかった場合
        if (!gapFound) {
            return false;
        }

        // 3. ノイズ除去判定
        // ★修正: 信号本体の幅が 200ns (0.2µs) 以上あるかチェック
        // (これで左右の空白2µsを足したトータル幅は 4.2µs 以上になる)
        double coreWidthNs = (double)(signalRightBin - signalLeftBin) * BIN_WIDTH_NS;

        if (coreWidthNs > minCoreSignalWidthNs)
        {
            // 合格: 積分範囲決定 (右端は +2µs)
            int integRightBin = signalRightBin + silenceThreshold;
            if (integRightBin > nBins) integRightBin = nBins;

            double sumWeight = 0;
            double sumCount = 0;
            double totalCounts = 0;

            // 積分実行 (左の空白 integLeftBin 〜 右の空白 integRightBin)
            for (int k = integLeftBin; k <= integRightBin; ++k)
            {
                double c = h->GetBinContent(k);
                double x = h->GetXaxis()->GetBinCenter(k);
                if (c > 0) {
                    sumWeight += (c * x);
                    sumCount += c;
                }
                totalCounts += c;
            }
            outCentroid = (sumCount > 0) ? (sumWeight / sumCount) : 0.0;
            outCount = totalCounts;
            return true;
        }
        else
        {
            // 不合格(ノイズ): この山をスキップして次へ
            searchHead = integLeftBin - 1;
            continue;
        }
    }
    return false;
}

// ----------------------------------------------------------------------------
// findBestShift: テンプレートマッチング用 (ライブラリ機能として温存)
// ----------------------------------------------------------------------------
int DetectorAnalyzer::findBestShift(TH1F *hTarget, TH1F *hTemplate, int binMin, int binMax)
{
    double minRes = 1.79769e+308;
    int bestS = SHIFT_CALC_ERROR;
    bool found = false;
    const int SEARCH_RANGE = 25; // 探索範囲

    double intTarget = hTarget->Integral(binMin, binMax);
    double intTemplate = hTemplate->Integral(binMin, binMax);

    if (intTarget <= 0 || intTemplate <= 0)
        return SHIFT_CALC_ERROR;

    double scale = intTemplate / intTarget;

    for (int s = -SEARCH_RANGE; s <= SEARCH_RANGE; ++s)
    {
        double res = calculateResidual(hTarget, hTemplate, s, binMin, binMax, scale);
        if (res >= 0 && res < minRes)
        {
            minRes = res;
            bestS = s;
            found = true;
        }
    }

    if (!found || std::abs(bestS) == SEARCH_RANGE)
        return SHIFT_CALC_ERROR;
    return bestS;
}

// ----------------------------------------------------------------------------
// calculateResidual: マッチング残差計算 (ライブラリ機能として温存)
// ----------------------------------------------------------------------------
double DetectorAnalyzer::calculateResidual(TH1F *hTarget, TH1F *hTemplate, int shiftBins, int binMin, int binMax, double scale)
{
    int tBinMin = binMin + shiftBins;
    int tBinMax = binMax + shiftBins;

    if (tBinMin < 1 || tBinMax > hTarget->GetNbinsX())
        return -1.0;

    double residualSum = 0.0;
    for (int i = binMin; i <= binMax; ++i)
    {
        int targetBin = i + shiftBins;
        double diff = (hTarget->GetBinContent(targetBin) * scale) - hTemplate->GetBinContent(i);
        residualSum += (diff * diff);
    }
    return residualSum;
}


// ----------------------------------------------------------------------------
// generatePDFReport: 修正版
//  - 重ね書きグラフ(緑)を「線と点 (PL)」で描画するように変更
//  - これにより視認性を向上 (飛んだ箇所も線で繋がるが許容する)
// ----------------------------------------------------------------------------
void DetectorAnalyzer::generatePDFReport()
{
    if (hRawToTMap.empty())
    {
        printLog("[WARN] No data available for PDF report.");
        return;
    }

    printLog("Generating PDF validation report (128 Strips + Summaries)...");
    std::string pdfName = "AnalysisReport.pdf";
    TCanvas *c1 = new TCanvas("c1", "Final Report", 1200, 800);
    c1->Print((pdfName + "(").c_str());

    std::string detNames[] = {"X1", "Y1", "X2", "Y2"};

    // ソート処理
    struct SortEntry
    {
        int type;
        int strip;
        int mod;
        int ch;
    };
    std::vector<SortEntry> sortedList;

    for (auto const &[sortKey, rawHist] : hRawToTMap)
    {
        if (!rawHist)
            continue;
        auto &cfg = detConfigMap_[sortKey];
        sortedList.push_back({cfg.detTypeID, cfg.strip, sortKey.first, sortKey.second});
    }
    std::sort(sortedList.begin(), sortedList.end(), [](const SortEntry &a, const SortEntry &b)
              {
        if (a.type != b.type) return a.type < b.type;
        return a.strip < b.strip; });

    c1->SetLogy(1);
    c1->SetRightMargin(0.12);

    for (const auto &entry : sortedList)
    {
        std::pair<int, int> plotKey = {entry.mod, entry.ch};

        c1->Clear();

        // 1. ヒストグラム描画
        TH1F *h = hRawToTMap[plotKey];
        h->SetTitle(Form("%s Strip%d (Mod%d Ch%d);ToT [ns];Counts",
                         detNames[entry.type].c_str(), entry.strip, entry.mod, entry.ch));
        h->GetXaxis()->SetRangeUser(0, MONITOR_HIST_MAX_TOT);
       h->SetLineColor(kBlue + 1);
        h->Draw();

        // ★追加: 補正後ヒストグラムの重ね書き
        int gIdx = (entry.type * 32) + (entry.strip - 1);
        if (hCorrectedToT_Global[gIdx])
        {
            TH1F* hCorr = hCorrectedToT_Global[gIdx];
            hCorr->SetLineColor(kRed);
            hCorr->SetLineWidth(1);
            hCorr->SetFillColorAlpha(kRed, 0.3); // 半透明の赤
            hCorr->Draw("same hist");
        }

        double yMax = h->GetMaximum() * 1.5;
        double yMin = 0.5;
        h->GetYaxis()->SetRangeUser(yMin, yMax);

        // 2. 赤枠
        int rMin = rangeMinMap_[plotKey];
        int rMax = rangeMaxMap_[plotKey];
        if (rMin > 0 && rMax > rMin)
        {
            double xMin = (double)rMin * BIN_WIDTH_NS;
            double xMax = (double)rMax * BIN_WIDTH_NS;
            TBox *box = new TBox(xMin, yMin, xMax, yMax);
            box->SetFillStyle(0);
            box->SetLineColor(kRed);
            box->SetLineWidth(2);
            box->Draw("same");
        }

        // 3. 重心推移 (緑の線 + 点)
        if (gainOverlayGraphs_.count(plotKey))
        {
            TGraph *gOriginal = gainOverlayGraphs_[plotKey];

            if (gOriginal && gOriginal->GetN() > 1)
            {
                TGraph *gVis = new TGraph();

                double x_dummy, tStart, tEnd;
                gOriginal->GetPoint(0, x_dummy, tStart);
                gOriginal->GetPoint(gOriginal->GetN() - 1, x_dummy, tEnd);

                double duration = tEnd - tStart;
                if (duration <= 0)
                    duration = 1.0;

                double logMin = std::log10(yMin);
                double logMax = std::log10(yMax);
                double logMid = (logMin + logMax) / 2.0;

                for (int i = 0; i < gOriginal->GetN(); ++i)
                {
                    double tot_ns, time_val;
                    gOriginal->GetPoint(i, tot_ns, time_val);

                    double ratio = (time_val - tStart) / duration;
                    double logY = logMid + (logMax - logMid) * ratio;
                    double plotY = std::pow(10.0, logY);

                    gVis->SetPoint(i, tot_ns, plotY);
                }

                gVis->SetMarkerStyle(20);
                gVis->SetMarkerSize(0.5);
                gVis->SetMarkerColor(kGreen + 2);
                gVis->SetLineColor(kGreen + 1); // 線の色設定
                gVis->SetLineWidth(1);          // 線の太さ設定

                // ★ここを修正: "PL same" (Points + Lines) に変更
                gVis->Draw("PL same");

                // 時間ラベル
                TLatex latex;
                latex.SetNDC(false);
                latex.SetTextSize(0.03);
                latex.SetTextColor(kGreen + 2);
                latex.SetTextAlign(12);

                double xRight = MONITOR_HIST_MAX_TOT * 1.02;
                double yMid = std::pow(10.0, logMid);

                latex.DrawLatex(xRight, yMid, "Time: 0 min");
                latex.DrawLatex(xRight, yMax, Form("Time: %.1f min", duration));

                TArrow *ar = new TArrow(xRight * 0.99, yMid, xRight * 0.99, yMax, 0.02, "|>");
                ar->SetLineColor(kGreen + 2);
                ar->Draw();
            }
        }

        c1->Print(pdfName.c_str());
    }

    // --- サマリーページ ---
    c1->Clear();
    c1->SetLogy(1);
    if (hDeltaT_Neighbor)
    {
        hDeltaT_Neighbor->Draw();
        c1->Print(pdfName.c_str());
    }
    if (hDeltaT_n8)
    {
        hDeltaT_n8->Draw();
        c1->Print(pdfName.c_str());
    }
    c1->SetLogy(0);
    if (hGlobalPIndexMap)
    {
        hGlobalPIndexMap->Draw();
        c1->Print(pdfName.c_str());
    }

    // Gain Trend
    for (int t = 0; t < 4; ++t)
    {
        c1->Clear();
        c1->SetGrid();
        TMultiGraph *mgG = new TMultiGraph();
        mgG->SetTitle(Form("Gain Peak History: %s;Time [min];Peak [ToT]", detNames[t].c_str()));
        for (auto const &[trendKey, g] : gainEvolutionGraphs_)
        {
            if (detConfigMap_[trendKey].detTypeID == t && g->GetN() > 0)
            {
                TGraph *gc = (TGraph *)g->Clone();
                gc->SetLineColor(detConfigMap_[trendKey].strip % 9 + 1);
                mgG->Add(gc);
            }
        }
        if (mgG->GetListOfGraphs())
        {
            mgG->Draw("AL");
            c1->Print(pdfName.c_str());
        }
    }

    // Rate Trend
    for (int t = 0; t < 4; ++t)
    {
        c1->Clear();
        c1->SetGrid();
        TMultiGraph *mgR = new TMultiGraph();
        mgR->SetTitle(Form("Rate History (CPM): %s;Time [min];Rate [CPM]", detNames[t].c_str()));
        for (auto const &[trendKey, g] : rateEvolutionGraphs_)
        {
            if (detConfigMap_[trendKey].detTypeID == t && g->GetN() > 0)
            {
                TGraph *rc = (TGraph *)g->Clone();
                rc->SetLineColor(detConfigMap_[trendKey].strip % 9 + 1);
                mgR->Add(rc);
            }
        }
        if (mgR->GetListOfGraphs())
        {
            mgR->Draw("AL");
            c1->Print(pdfName.c_str());
        }
    }

    c1->Print((pdfName + ")").c_str());
    delete c1;
    printLog("Full PDF report generated successfully.");
}

// 他の補助関数も全て維持
void DetectorAnalyzer::printFileTail(const std::string &filePath, long long currentOffset)
{
    std::ifstream file(filePath, std::ios::binary);
    if (!file.is_open())
        return;
    long long startPos = currentOffset > 64 ? currentOffset - 64 : 0;
    file.seekg(startPos);
    std::vector<unsigned char> buffer(128);
    file.read(reinterpret_cast<char *>(buffer.data()), buffer.size());
    std::streamsize bytesRead = file.gcount();
    std::cout << "\n========================================================" << std::endl;
    std::cout << " [HEX DUMP] File: " << filePath << std::endl;
    std::cout << " [HEX DUMP] Offset: " << startPos << " - " << (startPos + bytesRead) << std::endl;
    std::cout << "========================================================" << std::endl;
    for (int i = 0; i < bytesRead; ++i)
    {
        if (i % 16 == 0)
            std::cout << "\n " << std::setw(8) << std::setfill('0') << std::hex << (startPos + i) << ": ";
        std::cout << std::setw(2) << std::setfill('0') << std::hex << (int)buffer[i] << " ";
        if ((startPos + i) == currentOffset - 1)
            std::cout << "< ";
        else
            std::cout << " ";
    }
    std::cout << "\n========================================================\n"
              << std::dec << std::endl;
}

// --- クラスのメンバ関数としての printLog ---
void DetectorAnalyzer::printLog(const std::string &msg)
{
    std::cout << "\r\033[K" << std::flush;
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm *now_tm = std::localtime(&now_c);
    std::cout << "[" << std::put_time(now_tm, "%Y-%m-%d %H:%M:%S") << "] [INFO] " << msg << std::endl;
}

void DetectorAnalyzer::setupTree()
{
    if (!outputFile_ || !outputFile_->IsOpen())
        return;
    outputFile_->cd();

    if (tree_) delete tree_;
    tree_ = new TTree("tree", "Detector Events");
    tree_->Branch("type", &b_type, "type/B");
    tree_->Branch("strip", &b_strip, "strip/B");
    tree_->Branch("tot", &b_tot, "tot/I");
    tree_->Branch("time", &b_time, "time/L");
    
    // ★追加: メインTreeへの補正後ToTブランチ
    tree_->Branch("corrected_tot", &b_corrected_tot, "corrected_tot/D");

    const char *detNames[] = {"X1", "Y1", "X2", "Y2"};

    // ヒストグラムのビン数を定数から計算 (100000.0 / 1.0 = 100000 bins)
    int nBins = static_cast<int>(MONITOR_HIST_MAX_TOT / BIN_WIDTH_NS);

    for (int type = 0; type < 4; ++type)
    {
        for (int strip = 1; strip <= 32; ++strip)
        {
            int gIdx = (type * 32) + (strip - 1);

            std::string name = Form("hRaw_%s_Strip%02d", detNames[type], strip);
            std::string title = Form("%s Strip %02d;ToT [ns];Counts", detNames[type], strip);

            if (hRawToT_Global[gIdx]) delete hRawToT_Global[gIdx];
            hRawToT_Global[gIdx] = new TH1F(name.c_str(), title.c_str(), nBins, 0, MONITOR_HIST_MAX_TOT);

            std::string gName = Form("GCheck_%s_Strip%02d", detNames[type], strip);
            if (hGainCheck_Global[gIdx]) delete hGainCheck_Global[gIdx];
            hGainCheck_Global[gIdx] = new TH1F(gName.c_str(), title.c_str(), nBins, 0, MONITOR_HIST_MAX_TOT);

            // ★追加: 補正確認用ヒストグラム (赤色表示用) の初期化
            std::string cName = Form("hCorr_%s_Strip%02d", detNames[type], strip);
            std::string cTitle = Form("%s Strip %02d (Corrected);ToT [ns];Counts", detNames[type], strip);
            
            if (hCorrectedToT_Global[gIdx]) delete hCorrectedToT_Global[gIdx];
            hCorrectedToT_Global[gIdx] = new TH1F(cName.c_str(), cTitle.c_str(), nBins, 0, MONITOR_HIST_MAX_TOT);
            hCorrectedToT_Global[gIdx]->SetLineColor(kRed); 
        }
    }

    hDeltaT_Neighbor = new TH1F("DeltaT_Nearest", "Delta T (Nearest);ns;Counts", 100000, 0, 100000);
    hDeltaT_n8 = new TH1F("DeltaT_n8", "Delta T (N to N+7);ns;Counts", 100000, 0, 100000);
    hGlobalPIndexMap = new TH1F("GlobalMap", "Hit Map;Global Index;Counts", 128, 0, 128);

    printLog("setupTree: Initialized 128 Histograms (Dynamic 1ns/bin resolution).");
}

// ----------------------------------------------------------------------------
// 2. loadConversionFactors: マッピングのリンク (生成はしない)
// ----------------------------------------------------------------------------
void DetectorAnalyzer::loadConversionFactorsFromCSV(const std::string &csvFileName, const std::string &detName, int moduleID)
{
    std::ifstream file(csvFileName);
    if (!file.is_open())
        return;

    std::string line;
    std::getline(file, line); // ヘッダースキップ

    activeModuleIDs_.insert(moduleID);

    // このCSVファイルが担当する検出器タイプIDを特定
    int typeID = -1;
    if (detName == "X1")
        typeID = ID_X1;
    else if (detName == "Y1")
        typeID = ID_Y1;
    else if (detName == "X2")
        typeID = ID_X2;
    else if (detName == "Y2")
        typeID = ID_Y2;

    if (typeID == -1)
    {
        printLog("[WARN] Unknown detector name in CSV load: " + detName);
        return;
    }

    while (std::getline(file, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;
        while (std::getline(ss, cell, ','))
            row.push_back(cell);

        if (row.size() < 2)
            continue;

        try
        {
            int sysCh = std::stoi(row[0]);
            // Y検出器のチャンネルオフセット補正
            if (detName.find('Y') != std::string::npos)
                sysCh += Y_CHANNEL_OFFSET;
            if (sysCh >= MAX_SYS_CH)
                continue;

            int localStrip = std::stoi(row[1]); // 1-32
            if (localStrip < 1 || localStrip > 32)
                continue;

            // LUT登録
            ch_LUT[moduleID][sysCh].isValid = true;
            ch_LUT[moduleID][sysCh].detTypeID = typeID;
            ch_LUT[moduleID][sysCh].strip = localStrip;

            // ★リンク処理: 先行生成されたGlobal配列から、該当するポインタを探してLUTに繋ぐ
            int gIdx = (typeID * 32) + (localStrip - 1);

            if (gIdx >= 0 && gIdx < 128)
            {
                // LUTにポインタをコピー（これでprocessChunkから高速アクセス可能）
                hRawToT_LUT[moduleID][sysCh] = hRawToT_Global[gIdx];
                hGainCheck_LUT[moduleID][sysCh] = hGainCheck_Global[gIdx];

                // Mapにも登録（保存やPDF生成で使う場合用）
                hRawToTMap[{moduleID, sysCh}] = hRawToT_Global[gIdx];
                gainCheckHists_[{moduleID, sysCh}] = hGainCheck_Global[gIdx];

                // コンフィグマップにも登録
                detConfigMap_[{moduleID, sysCh}] = {typeID, localStrip, detName, -1};
            }
        }
        catch (...)
        {
            continue;
        }
    }
    printLog("Mapped " + detName + " to Pre-allocated Histograms for Module " + std::to_string(moduleID));
}

// --- あかりさんのルールに基づくラン識別子抽出 ---
std::string DetectorAnalyzer::parseRunPrefix(const std::string &fullPath)
{
    std::string fileName = std::filesystem::path(fullPath).filename().string();
    size_t firstUnderscore = fileName.find('_');
    if (firstUnderscore != std::string::npos)
    {
        return fileName.substr(0, firstUnderscore); // "PSD000001"
    }
    return fileName;
}

// --- フォルダ + PSD でランを特定する署名 ---
std::string DetectorAnalyzer::getRunSignature(const std::string &f)
{
    std::filesystem::path p(f);
    return p.parent_path().string() + "###" + parseRunPrefix(f);
}

double DetectorAnalyzer::calculatePeakAreaCovell(TH1F *hist, double peakPos)
{
    if (!hist || peakPos <= 0)
        return 0.0;
    int binCenter = hist->FindBin(peakPos);
    int halfWidth = COVELL_WINDOW_NS / 20;
    int binMin = binCenter - halfWidth;
    int binMax = binCenter + halfWidth;
    if (binMin < 1)
        binMin = 1;
    if (binMax > hist->GetNbinsX())
        binMax = hist->GetNbinsX();
    if (binMax <= binMin)
        return 0.0;
    double grossArea = hist->Integral(binMin, binMax);
    double countLeft = hist->GetBinContent(binMin);
    double countRight = hist->GetBinContent(binMax);
    double numBins = (double)(binMax - binMin + 1);
    double bgArea = (countLeft + countRight) * numBins / 2.0;
    double netArea = grossArea - bgArea;
    return (netArea < 0) ? 0.0 : netArea;
}

unsigned long long DetectorAnalyzer::getSafeTime(const std::map<int, unsigned long long> &lastTimes)
{
    unsigned long long minT = std::numeric_limits<unsigned long long>::max();
    if (lastTimes.empty())
        return 0;
    for (auto const &[id, t] : lastTimes)
    {
        if (t < minT)
            minT = t;
    }
    return minT;
}

short DetectorAnalyzer::getRunID(const std::string &prefix)
{
    if (runPrefixToId_.count(prefix))
        return runPrefixToId_[prefix];
    short newId = (short)runIdToPrefix_.size();
    runIdToPrefix_.push_back(prefix);
    runPrefixToId_[prefix] = newId;
    return newId;
}

void DetectorAnalyzer::loadAndSortFiles(const std::string &directory, std::map<int, std::deque<std::string>> &fileQueues)
{
    // 1. ファイル走査
    struct FileInfo
    {
        int runID;
        int modID;
        int fileSeq;
        std::string path;
        std::string name;
    };
    std::vector<FileInfo> tempAllFiles;
    totalDataSize_ = 0;

    if (!fs::exists(directory))
    {
        printLog("[ERROR] Directory not found: " + directory);
        return;
    }

    auto parseAndAdd = [&](const fs::path &p)
    {
        if (p.extension() != ".edb")
            return;
        std::string fullPath = p.string();
        if (fullPath.find("original") != std::string::npos)
            return;

        std::string fName = p.stem().string();
        int fRun = 0, fMod = 0, fSeq = 0;
        if (std::sscanf(fName.c_str(), "PSD%d_00_%d_%d", &fRun, &fMod, &fSeq) == 3)
        {
            tempAllFiles.push_back({fRun, fMod, fSeq, fullPath, fName});
            totalDataSize_ += fs::file_size(p);
        }
    };

    printLog("Scanning directory: " + directory);
    for (const auto &entry : fs::directory_iterator(directory))
    {
        if (entry.is_directory())
        {
            std::string dName = entry.path().filename().string();
            if (dName == "original")
                continue;
            for (const auto &sub : fs::directory_iterator(entry.path()))
                parseAndAdd(sub.path());
        }
        else
        {
            parseAndAdd(entry.path());
        }
    }

    // 2. ソート (RunID -> Seq -> Mod)
    std::sort(tempAllFiles.begin(), tempAllFiles.end(), [](const FileInfo &a, const FileInfo &b)
              {
        if (a.runID != b.runID) return a.runID < b.runID;
        if (a.fileSeq != b.fileSeq) return a.fileSeq < b.fileSeq;
        return a.modID < b.modID; });

    // 3. キューへ投入
    for (const auto &info : tempAllFiles)
    {
        fileQueues[info.modID].push_back(info.path);
    }

    // ========================================================================
    //  左右並列リスト表示 (警告なし・シンプル版)
    // ========================================================================
    std::map<int, std::map<int, std::map<int, std::string>>> visualMap;
    for (const auto &info : tempAllFiles)
    {
        visualMap[info.runID][info.fileSeq][info.modID] = info.name;
    }

    printLog("========================================================================================");
    printLog(" [INFO] File Queue Content (Side-by-Side)");
    printLog("========================================================================================");
    std::cout << " RunID  | Seq | Module 0                       | Module 1                       " << std::endl;
    std::cout << "--------+-----+--------------------------------+--------------------------------" << std::endl;

    for (const auto &[runID, seqMap] : visualMap)
    {
        for (const auto &[seq, modMap] : seqMap)
        {
            std::string name0 = (modMap.count(0)) ? modMap.at(0) : "";
            std::string name1 = (modMap.count(1)) ? modMap.at(1) : "";

            // 長すぎる場合は省略
            if (name0.length() > 30)
                name0 = ".." + name0.substr(name0.length() - 28);
            if (name1.length() > 30)
                name1 = ".." + name1.substr(name1.length() - 28);

            std::cout << " " << std::setw(6) << runID
                      << " | " << std::setw(3) << seq
                      << " | " << std::left << std::setw(30) << name0
                      << " | " << std::left << std::setw(30) << name1 << std::endl;
        }
        std::cout << "--------+-----+--------------------------------+--------------------------------" << std::endl;
    }

    printLog("Queue Setup Complete.");
    printLog("Mod 0 Files: " + std::to_string(fileQueues[0].size()));
    printLog("Mod 1 Files: " + std::to_string(fileQueues[1].size()));
    printLog("========================================================================================");
}

void DetectorAnalyzer::loadOfflineStripList(const std::string &csvFileName)
{
    std::ifstream file(csvFileName);
    if (!file.is_open())
        return;
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;
        while (std::getline(ss, cell, ','))
            row.push_back(cell);
        if (row.size() < 3)
            continue;
        try
        {
            int detType = std::stoi(row[0]);
            int stripNo = std::stoi(row[1]);
            int flag = std::stoi(row[2]);
            if (flag == 1)
                offlineStrips_.insert({detType, stripNo});
        }
        catch (...)
        {
        }
    }
}

void DetectorAnalyzer::setTimeLimitMinutes(double minutes)
{
    if (minutes <= 0)
    {
        analysisDuration_ns_ = std::numeric_limits<unsigned long long>::max();
        printLog(">>> Time limit set to: UNLIMITED");
    }
    else
    {
        analysisDuration_ns_ = static_cast<unsigned long long>(minutes * 60.0 * 1.0e9);
        std::stringstream ss;
        ss << ">>> Time limit set to: " << std::fixed << std::setprecision(2) << minutes << " minutes";
        printLog(ss.str());
    }
}

// ----------------------------------------------------------------------------
// fillFluxConserving: 1nsのビン幅を保ったまま、スケール変換後のカウントを分配する
// ----------------------------------------------------------------------------
void DetectorAnalyzer::fillFluxConserving(TH1F* h, double rawToT, double scale)
{
    // 元の区間 [rawToT, rawToT + 1.0) を scale 倍する
    double low = rawToT * scale;
    double high = (rawToT + 1.0) * scale;
    
    // 変換後の幅 (線形変換なので常に scale と等しい)
    double width = high - low; 
    if (width <= 0.0) return;

    // 影響するビンの範囲
    int binStart = (int)std::floor(low);
    int binEnd = (int)std::floor(high);

    // 各ビンへの分配
    for (int b = binStart; b <= binEnd; ++b)
    {
        double binLo = (double)b;
        double binHi = (double)(b + 1);

        double overlapLo = std::max(low, binLo);
        double overlapHi = std::min(high, binHi);
        double overlapLen = overlapHi - overlapLo;

        if (overlapLen > 0.0)
        {
            double weight = overlapLen / width;
            h->Fill(binLo + 0.5, weight);
        }
    }
}