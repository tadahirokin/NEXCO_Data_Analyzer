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

namespace fs = std::filesystem;

// ----------------------------------------------------------------------------
// グローバル変数 & ログ関数
// ----------------------------------------------------------------------------
static DetectorAnalyzer* g_currentAnalyzer = nullptr;

void printLog(const std::string& msg) {
    std::cout << "\r\033[K" << std::flush;
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm* now_tm = std::localtime(&now_c);
    std::cout << "[" << std::put_time(now_tm, "%Y-%m-%d %H:%M:%S") << "] [INFO] " << msg << std::endl;

    if (g_currentAnalyzer && g_currentAnalyzer->isAnalysisStarted_) {
        g_currentAnalyzer->printSearchStatus(); 
    }
}

// ----------------------------------------------------------------------------
// コンストラクタ
// ----------------------------------------------------------------------------
DetectorAnalyzer::DetectorAnalyzer(int timeWindow_ns, const std::string& outputFileName)
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
      abortCurrentRun_(false), // 新規フラグ初期化
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
      b_i_type(nullptr), b_i_strip(nullptr), b_i_tot(nullptr), b_i_time(nullptr),
      globalEventID_Ideal_(0)
{
    g_currentAnalyzer = this;
    gROOT->SetBatch(kTRUE); 

    analysisStartTime_ = std::chrono::system_clock::now();

    if (!outputFile_ || outputFile_->IsZombie()) {
        throw std::runtime_error("[ERROR] Failed to create output ROOT file.");
    }
    
    for(int m=0; m<MAX_MODULES; ++m) {
        for(int s=0; s<MAX_SYS_CH; ++s) {
            ch_LUT[m][s].isValid = false;
            hRawToT_LUT[m][s] = nullptr;
            hGainCheck_LUT[m][s] = nullptr;
        }
    }

    setupTree(); 

    struct stat st = {}; 
    if (stat("./gainshift", &st) == -1) {
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

void DetectorAnalyzer::setMuonAnalysisMode(bool enable) {
    enableMuonAnalysis_ = enable;
    if (enable) {
        printLog("[MODE] Muon Analysis ENABLED (Clustering Active, Slower)");
    } else {
        printLog("[MODE] Fast Gain Check ONLY (Muon Extraction Disabled)");
    }
}

DetectorAnalyzer::~DetectorAnalyzer() {
    printLog("Shutting down analyzer...");

    // 1. 最後のゲイン解析を実行
    if (currentEventTime_ns_ > lastGlobalAnalysisTime_ns_) {
        analyzeGainShift(currentEventTime_ns_);
    }

    // 2. 有効時間の最終計算とPDFレポートの生成
    calculateEffectiveTime();
    generatePDFReport(); 

    // 3. ROOTファイルの保存
    if (outputFile_ && outputFile_->IsOpen()) {
        outputFile_->cd();
        if (tree_) tree_->Write();
        
        // 解析結果（ヒストグラム）を保存
        auto dirHist = outputFile_->mkdir("Histograms"); 
        if (dirHist) {
            dirHist->cd();
            for (auto const& [key, h] : hRawToTMap) if (h) h->Write();
        }

        outputFile_->cd();
        auto dG = outputFile_->mkdir("GainHistory");
        if (dG) {
            dG->cd();
            for (auto const& [key, g] : gainEvolutionGraphs_) if (g && g->GetN() > 0) g->Write();
        }

        outputFile_->Close();
        delete outputFile_;
        outputFile_ = nullptr;
    }

    // 4. ミュオン抽出結果の保存
    if (idealFile_ && idealFile_->IsOpen()) {
        idealFile_->cd();
        if (idealTree_) idealTree_->Write();
        idealFile_->Close();
        delete idealFile_;
    }

    if (gainLogCsv_.is_open()) gainLogCsv_.close();
    if (rateLogCsv_.is_open()) rateLogCsv_.close();
}


std::pair<unsigned long long, bool> DetectorAnalyzer::readEventsFromFile(const std::string& fileName, int modID, std::vector<Event>& rawEvents, long long& offset) {
    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs.is_open()) return {lastT0_[modID], true};

    ifs.seekg(0, std::ios::end);
    long long fileSize = static_cast<long long>(ifs.tellg());
    if (offset >= fileSize) return {lastT0_[modID], true};

    const size_t CHUNK_SIZE = 64 * 1024 * 1024;

    // --- Phase 1: 物理アライメント (10パケット検証) ---
    if (!hasDataAligned_[modID]) {
        ifs.seekg(offset, std::ios::beg);
        std::vector<char> scanBuf(CHUNK_SIZE);
        ifs.read(scanBuf.data(), CHUNK_SIZE);
        size_t count = static_cast<size_t>(ifs.gcount());
        unsigned char* buf = reinterpret_cast<unsigned char*>(scanBuf.data());

        for (size_t i = 0; i + 80 < count; ++i) {
            if (buf[i] == 0x69 || buf[i] == 0x6a) {
                bool match = true;
                for (int k = 1; k < 10; ++k) {
                    if (buf[i + k * 8] != 0x69 && buf[i + k * 8] != 0x6a) { match = false; break; }
                }
                if (match) {
                    hasDataAligned_[modID] = true;
                    offset += static_cast<long long>(i);
                    break;
                }
            }
        }
        if (!hasDataAligned_[modID]) return {lastT0_[modID], (count < CHUNK_SIZE)};
    }

    // --- Phase 2: ストリーミング読み込み ---
    ifs.seekg(offset, std::ios::beg);
    std::vector<char> buffer(CHUNK_SIZE);
    ifs.read(buffer.data(), CHUNK_SIZE);
    size_t readCount = static_cast<size_t>(ifs.gcount());
    bool isEOF = (offset + static_cast<long long>(readCount) >= fileSize);
    unsigned char* buf = reinterpret_cast<unsigned char*>(buffer.data());

    rawEvents.reserve(rawEvents.size() + (readCount / 8));

    unsigned long long currentBaseTime = lastT0_[modID];
    size_t i = 0;
    size_t lastT0_offset_in_buf = 0;
    size_t valid_event_count_at_lastT0 = rawEvents.size();
    bool foundAnyT0 = false;

    while (i + 8 <= readCount) {
        if (buf[i] == 0x69) {
            unsigned char* p = &buf[i + 1];
            unsigned long long t = ((unsigned long long)p[0] << 24 | (unsigned long long)p[1] << 16 | 
                                   (unsigned long long)p[2] << 8 | (unsigned long long)p[3]) * 1000000000ULL + 
                                   ((unsigned long long)p[4] << 2 | (unsigned long long)(p[5] & 0xC0) >> 6) * 1000000ULL + 
                                   ((unsigned long long)(p[5] & 0x3F) << 2 | (unsigned long long)(p[6] & 0xC0) >> 6);

            if (t <= lastT0_[modID] && lastT0_[modID] != 0) { i += 8; continue; }

            currentBaseTime = t;
            foundAnyT0 = true;
            lastT0_offset_in_buf = i;
            valid_event_count_at_lastT0 = rawEvents.size();
            i += 8;
        } 
        else if (buf[i] == 0x6a) {
            if (currentBaseTime > 0) {
                unsigned char* p = &buf[i + 1];
                unsigned int tof = ((unsigned int)p[0] << 16 | (unsigned int)p[1] << 8 | (unsigned int)p[2]);
                unsigned int pw  = ((unsigned int)p[3] << 12 | (unsigned int)p[4] << 4 | (unsigned int)(p[5] & 0xF0) >> 4);
                int sys = static_cast<int>(((p[5] & 0x0F) << 8) | p[6]);
                unsigned long long evTime = currentBaseTime + static_cast<unsigned long long>(tof) - static_cast<unsigned long long>(pw);
                rawEvents.push_back({0, 0, modID, evTime, static_cast<int>(pw), sys});
            }
            i += 8;
        } 
        else { hasDataAligned_[modID] = false; break; }
    }

    // --- Phase 3: TOアンカーによる処理範囲の確定 (EOF時も共通) ---
    if (foundAnyT0) {
        // 最後のT0以降のデータは常に切り捨てる (次のファイルで重複して読むため)
        rawEvents.resize(valid_event_count_at_lastT0);
        
        if (isEOF) {
            offset = fileSize; // ファイルを閉じるために最後まで進める
        } else {
            offset += static_cast<long long>(lastT0_offset_in_buf); // 次回はこのT0から再開
        }
        return {currentBaseTime, isEOF};
    } else {
        // T0が見つからなかった場合
        offset += static_cast<long long>(readCount);
        return {currentBaseTime, isEOF};
    }
}


void DetectorAnalyzer::printSearchStatus() {
    double liveTimeMin = (double)totalEffectiveTimeNs_ / 1.0e9 / 60.0;
    double progress = (analysisDuration_ns_ > 0 && analysisDuration_ns_ != ULLONG_MAX) 
                      ? (double)totalEffectiveTimeNs_ / (double)analysisDuration_ns_ * 100.0 : 0.0;
    std::cout << "\r\033[K" << "Progress: [" << std::fixed << std::setprecision(1) << progress << "%] "
              << "Live Time: " << liveTimeMin << " min" << std::flush;
}

void DetectorAnalyzer::processBinaryFiles(std::map<int, std::deque<string>>& fileQueues) {
    const size_t PER_MOD_CAP = 10000000;
    std::map<int, std::deque<Event>> moduleBuffers;
    std::map<int, unsigned long long> moduleLastTime;
    unsigned long long lastProcessedTime = 0;
    auto lastUIDraw = std::chrono::system_clock::now();

    activeModuleIDs_.clear();
    for(auto const& [m, q] : fileQueues) if(!q.empty()) activeModuleIDs_.insert(m);

    for (int m : activeModuleIDs_) {
        lastT0_[m] = 0; hasFoundT0_[m] = false; hasDataAligned_[m] = false;
        moduleLastTime[m] = 0; currentFileOffsets_[m] = 0;
    }

    string currentRunSignature = "";

    while (true) {
        string nextFileSig = "";
        bool allQueuesEmpty = true;
        for (auto const& [m, q] : fileQueues) {
            if (!q.empty()) {
                allQueuesEmpty = false;
                if (nextFileSig.empty()) nextFileSig = getRunSignature(q.front());
            }
        }
        if (allQueuesEmpty && moduleBuffers.empty()) break;

        // --- ランの切り替え ---
        if (currentRunSignature != nextFileSig && !nextFileSig.empty()) {
            if (outputFile_) { outputFile_->Write(); outputFile_->Close(); delete outputFile_; outputFile_ = nullptr; }
            currentRunSignature = nextFileSig;
            string flat = currentRunSignature; for (char &c : flat) if (c == '/' || c == ':' || c == '\\') c = '_';
            outputFile_ = new TFile((flat + ".root").c_str(), "RECREATE");
            setupTree();
            
            lastProcessedTime = 0; // 進捗の起点をリセット
            for (int m : activeModuleIDs_) {
                hasDataAligned_[m] = false; hasFoundT0_[m] = false;
                lastT0_[m] = 0; currentFileOffsets_[m] = 0;
                moduleBuffers[m].clear(); moduleLastTime[m] = 0;
            }
        }

        // --- 補充ロジック ---
        int lagMod = -1; unsigned long long minT = ULLONG_MAX;
        for (int m : activeModuleIDs_) {
            if (moduleLastTime[m] < minT) { minT = moduleLastTime[m]; lagMod = m; }
        }

        if (lagMod != -1 && moduleBuffers[lagMod].size() < PER_MOD_CAP) {
            std::vector<Event> temp;
            auto res = readEventsFromFile(fileQueues[lagMod].front(), lagMod, temp, currentFileOffsets_[lagMod]);
            if (!temp.empty()) {
                moduleBuffers[lagMod].insert(moduleBuffers[lagMod].end(), 
                                             std::make_move_iterator(temp.begin()), 
                                             std::make_move_iterator(temp.end()));
            }
            moduleLastTime[lagMod] = res.first;
            if (res.second) { fileQueues[lagMod].pop_front(); currentFileOffsets_[lagMod] = 0; }
        }

        // --- SafeTime 同期と Progress 計算の修正 ---
        unsigned long long safeTime = getSafeTime(moduleLastTime);

        if (safeTime > lastProcessedTime) {
            // 初動の同期: 0 から実時刻へジャンプして起点を確定
            if (lastProcessedTime == 0) {
                lastProcessedTime = safeTime;
            } else {
                // 確定した時間幅を実効解析時間に加算
                totalEffectiveTimeNs_ += (safeTime - lastProcessedTime);
                lastProcessedTime = safeTime;
            }
            currentEventTime_ns_ = safeTime;

            // --- 解析フェーズ ---
            std::vector<Event> chunk;
            for (int m : activeModuleIDs_) {
                while (!moduleBuffers[m].empty() && moduleBuffers[m].front().eventTime_ns <= safeTime) {
                    chunk.push_back(std::move(moduleBuffers[m].front()));
                    moduleBuffers[m].pop_front();
                }
            }

            if (!chunk.empty()) {
                std::sort(chunk.begin(), chunk.end(), [](const Event& a, const Event& b) {
                    return a.eventTime_ns < b.eventTime_ns;
                });

                for (auto& ev : chunk) {
                    if (ev.module < MAX_MODULES && ev.sysCh < MAX_SYS_CH && ch_LUT[ev.module][ev.sysCh].isValid) {
                        int detType = ch_LUT[ev.module][ev.sysCh].detTypeID;
                        int strip = ch_LUT[ev.module][ev.sysCh].strip;
                        if (hRawToTMap.count({detType, strip})) hRawToTMap[{detType, strip}]->Fill(static_cast<double>(ev.tot));
                        
                        b_type = static_cast<Char_t>(detType);
                        b_strip = static_cast<Char_t>(strip);
                        b_tot = static_cast<Int_t>(ev.tot);
                        b_time = static_cast<Long64_t>(ev.eventTime_ns);
                        tree_->Fill();
                    }
                }
            }
        }

        // --- UI表示更新 ---
        if (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - lastUIDraw).count() > 500) {
            printSearchStatus();
            lastUIDraw = std::chrono::system_clock::now();
        }

        if (totalEffectiveTimeNs_ >= analysisDuration_ns_) break;
    }
}


unsigned long long DetectorAnalyzer::findNextT0(const std::string& fileName, int modID, long long& offset) {
    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs.is_open()) return 0;

    ifs.seekg(0, std::ios::end);
    long long fileSize = ifs.tellg();
    
    if (offset >= fileSize) return 0;

    // --- Phase 1: アライメント探索 (1バイトずつ走査) ---
    if (!hasDataAligned_[modID]) {
        ifs.seekg(offset, std::ios::beg);
        const size_t SCAN_BUF_SIZE = 64 * 1024 * 1024; // 64MB
        std::vector<char> scanBuf(SCAN_BUF_SIZE);
        ifs.read(scanBuf.data(), SCAN_BUF_SIZE);
        size_t count = ifs.gcount();
        unsigned char* buf = (unsigned char*)scanBuf.data();

        for (size_t i = 0; i + 80 < count; ++i) {
            unsigned char h = buf[i];
            if (h == 0x69 || h == 0x6a) {
                // 10パケット検証
                bool match = true;
                for (int k = 1; k < 10; ++k) {
                    if (buf[i + k*8] != 0x69 && buf[i + k*8] != 0x6a) { match = false; break; }
                }
                if (match) {
                    hasDataAligned_[modID] = true;
                    offset += i; // 1バイトずつの走査で確定した位置
                    break;
                }
            }
        }
    }

    if (!hasDataAligned_[modID]) return 0;

    // --- Phase 2: T0パケット探索 (確定した枠から8バイト読み) ---
    ifs.seekg(offset, std::ios::beg);
    unsigned char packet[8];
    while (ifs.read((char*)packet, 8)) {
        if (packet[0] == 0x69) {
            unsigned char* p = &packet[1];
            unsigned long long t = ((unsigned long long)p[0]<<24|p[1]<<16|p[2]<<8|p[3])*1000000000ULL + 
                                   ((unsigned long long)p[4]<<2|(p[5]&0xC0)>>6)*1000000ULL + 
                                   ((p[5]&0x3F)<<2|(p[6]&0xC0)>>6);
            // offsetをこのT0パケットの位置に更新（これを親が保持する）
            return t;
        }
        offset += 8;
    }

    return 0;
}

bool DetectorAnalyzer::processChunk(const std::vector<Event>& sortedEvents) {
    if (sortedEvents.empty()) return true;

    for (size_t i = 0; i < sortedEvents.size(); ++i) {
        const auto& e = sortedEvents[i];
        
        // 1. 死活監視 (5秒ルール)
        if (ignoredMonitorChannels_.count({e.module, e.sysCh})) {
            channelLastRefHitTime_[e.module] = e.eventTime_ns;
        }
        if (channelLastRefHitTime_[e.module] > 0) {
            if (e.eventTime_ns > channelLastRefHitTime_[e.module] + MODULE_DEATH_TIMEOUT_NS) {
                isCurrentRunDead_ = true;
                return false; 
            }
        }

        // 2. 有効時間の積算 (絶対時刻ベース)
        currentEventTime_ns_ = e.eventTime_ns;
        if (prevEventTimeForEff_ > 0) {
            unsigned long long dt_ns = e.eventTime_ns - prevEventTimeForEff_;
            if (dt_ns < 1000000000ULL) totalEffectiveTimeNs_ += dt_ns; 
        }
        prevEventTimeForEff_ = e.eventTime_ns;

        // 3. Treeへの書き込み (絶対時刻)
        b_type = static_cast<Char_t>(e.type);
        b_strip = static_cast<Char_t>(e.strip);
        b_tot = static_cast<Int_t>(e.tot);
        b_time = static_cast<Long64_t>(e.eventTime_ns); // globalStartを引かない
        if (tree_) tree_->Fill();

        // 4. 各種統計ヒストグラムのFill
        if (hGainCheck_LUT[e.module][e.sysCh]) hGainCheck_LUT[e.module][e.sysCh]->Fill((double)e.tot);
        if (hRawToT_LUT[e.module][e.sysCh]) hRawToT_LUT[e.module][e.sysCh]->Fill((double)e.tot);
        
        int pIdx = ch_LUT[e.module][e.sysCh].pIndex;
        if (hGlobalPIndexMap && pIdx >= 0) hGlobalPIndexMap->Fill((double)pIdx);

        // 隣接イベント間隔 (DeltaT Neighbor)
        if (hDeltaT_Neighbor && i > 0) {
            hDeltaT_Neighbor->Fill((double)(e.eventTime_ns - sortedEvents[i-1].eventTime_ns));
        }
        // 8イベント間隔
        if (i + 7 < sortedEvents.size() && hDeltaT_n8) {
            double dt8 = (double)(sortedEvents[i+7].eventTime_ns - e.eventTime_ns);
            if (dt8 >= 0 && dt8 < 1000000) hDeltaT_n8->Fill(dt8);
        }

        // 5. ミュオン解析モード時のみ：バッファに蓄積
        if (enableMuonAnalysis_) {
            analysisBuffer_.push_back(e);
        }
    }

    // 6. ミュオン抽出処理の実行
    if (enableMuonAnalysis_) {
        processEventExtraction();
    }

    return true;
}

void DetectorAnalyzer::processEventExtraction() {
    if (analysisBuffer_.empty()) return;
    while (analysisBuffer_.size() >= 8) {
        double w = 0;
        bool isA = checkConditionA(analysisBuffer_, 0, w);
        bool isB = isA ? false : checkConditionB(analysisBuffer_, 0, w);

        if (isA || isB) {
            b_i_runID = currentProcessingRunID_;
            b_i_eventID = globalEventID_Ideal_++;
            b_i_nHits = isA ? 8 : 7;
            b_i_width = w;
            b_i_type->clear(); b_i_strip->clear(); b_i_tot->clear(); b_i_time->clear();
            for (int k = 0; k < b_i_nHits; ++k) {
                const auto& e = analysisBuffer_[k];
                b_i_type->push_back(ch_LUT[e.module][e.sysCh].detTypeID);
                b_i_strip->push_back(ch_LUT[e.module][e.sysCh].strip);
                b_i_tot->push_back(e.tot);
                b_i_time->push_back((double)e.eventTime_ns);
            }
            if (idealTree_) idealTree_->Fill();
            for (int k = 0; k < b_i_nHits; ++k) analysisBuffer_.pop_front();
        } else {
            analysisBuffer_.pop_front();
        }
    }
}

// 条件A: 全層で隣り合う2つを含む (8ヒット)
bool DetectorAnalyzer::checkConditionA(const std::deque<Event>& buf, size_t startIdx, double& outWidth) {
    if (startIdx + 8 > buf.size()) return false;
    outWidth = (double)(buf[startIdx+7].eventTime_ns - buf[startIdx].eventTime_ns);

    std::map<int, std::vector<int>> sMap;
    for(size_t i=0; i<8; ++i) {
        const auto& e = buf[startIdx+i];
        int type = ch_LUT[e.module][e.sysCh].detTypeID;
        int strip = ch_LUT[e.module][e.sysCh].strip;
        sMap[type].push_back(strip);
    }
    // X1, Y1, X2, Y2 全てで隣接ペアを持つか
    for(int type=0; type<4; ++type) {
        if (!hasAdjacentPair(sMap[type])) return false;
    }
    return true;
}

// 条件B: Y1欠損考慮 (7ヒット)
bool DetectorAnalyzer::checkConditionB(const std::deque<Event>& buf, size_t startIdx, double& outWidth) {
    if (startIdx + 7 > buf.size()) return false;
    outWidth = (double)(buf[startIdx+6].eventTime_ns - buf[startIdx].eventTime_ns);
    
    std::map<int, std::vector<int>> sMap;
    for(size_t i=0; i<7; ++i) {
        const auto& e = buf[startIdx+i];
        int type = ch_LUT[e.module][e.sysCh].detTypeID;
        int strip = ch_LUT[e.module][e.sysCh].strip;
        sMap[type].push_back(strip);
    }

    if (!hasAdjacentPair(sMap[ID_X1])) return false;
    if (!hasAdjacentPair(sMap[ID_X2])) return false;
    if (!hasAdjacentPair(sMap[ID_Y2])) return false;
    
    if (sMap[ID_Y1].size() != 1) return false;
    return isAdjacentToOffline(sMap[ID_Y1][0], ID_Y1);
}

bool DetectorAnalyzer::hasAdjacentPair(const std::vector<int>& strips) {
    if (strips.size() < 2) return false;
    std::vector<int> s = strips;
    std::sort(s.begin(), s.end());
    for(size_t i=0; i<s.size()-1; ++i) {
        if (s[i+1] - s[i] == 1) return true;
    }
    return false;
}

bool DetectorAnalyzer::isAdjacentToOffline(int strip, int detTypeID) {
    if (offlineStrips_.count({detTypeID, strip - 1})) return true;
    if (offlineStrips_.count({detTypeID, strip + 1})) return true;
    return false;
}

// ----------------------------------------------------------------------------
// ゲイン解析 (テンプレートマッチング)
// ----------------------------------------------------------------------------
bool DetectorAnalyzer::analyzeGainShift(unsigned long long currentTime_ns) {
    double timeMin = currentTime_ns / 60.0e9;

    if (!isCsvHeaderWritten_) {
        std::string detNames[] = {"X1", "Y1", "X2", "Y2"};
        gainLogCsv_ << "Time_ns"; rateLogCsv_ << "Time_ns";
        for (int type = 0; type < 4; ++type) {
            for (int s = 1; s <= 32; ++s) {
                std::string header = "," + detNames[type] + "_" + std::to_string(s);
                gainLogCsv_ << header; rateLogCsv_ << header;
            }
        }
        gainLogCsv_ << "\n"; rateLogCsv_ << "\n";
        isCsvHeaderWritten_ = true;
    }

    std::vector<double> currentStepPeaks(128, 0.0);
    std::vector<double> currentStepRates(128, 0.0);

    for (auto const& [key, hist] : gainCheckHists_) {
        auto& cfg = detConfigMap_[key];
        if (!hist) continue;

        // テンプレート作成
        if (templateHists_.find(key) == templateHists_.end()) {
            int peakBin = findRightMostPeak(hist);
            int rMin = 0, rMax = 0;
            if (peakBin > 0) determineIntegrationRange(hist, peakBin, rMin, rMax);

            if (rMin > 0 && rMax > rMin && hist->GetEntries() > 50) {
                rangeMinMap_[key] = rMin;
                rangeMaxMap_[key] = rMax;
                initialPeakPosMap_[key] = hist->GetBinCenter(peakBin);
                cumulativeShiftNsMap_[key] = 0.0;

                TH1F* tHist = (TH1F*)hist->Clone(Form("Template_Mod%d_Ch%d", key.first, key.second));
                tHist->SetDirectory(0);
                templateHists_[key] = tHist;

                if (!gainEvolutionGraphs_[key]) gainEvolutionGraphs_[key] = new TGraph();
                gainEvolutionGraphs_[key]->SetPoint(0, timeMin, initialPeakPosMap_[key]);
            }
        } 
        // テンプレートマッチング
        else if (rangeMinMap_[key] > 0) {
            TH1F* tHist = templateHists_[key];
            int deltaBin = findBestShift(hist, tHist, rangeMinMap_[key], rangeMaxMap_[key], cfg.cachedName);
            
            if (deltaBin != SHIFT_CALC_ERROR) {
                rangeMinMap_[key] += deltaBin;
                rangeMaxMap_[key] += deltaBin;
                cumulativeShiftNsMap_[key] += (deltaBin * 20.0);

                double absPeakPos = initialPeakPosMap_[key] + cumulativeShiftNsMap_[key];
                if (!gainEvolutionGraphs_[key]) gainEvolutionGraphs_[key] = new TGraph();
                gainEvolutionGraphs_[key]->SetPoint(gainEvolutionGraphs_[key]->GetN(), timeMin, absPeakPos);
                
                tHist->Reset();
                tHist->Add(hist);
            }
        }

        if (rangeMinMap_[key] > 0) {
            int idx = cfg.detTypeID * 32 + (cfg.strip - 1);
            if (idx >= 0 && idx < 128) {
                currentStepPeaks[idx] = initialPeakPosMap_[key] + cumulativeShiftNsMap_[key];
                double counts = hist->Integral(rangeMinMap_[key], rangeMaxMap_[key]);
                currentStepRates[idx] = counts / (600.0) * 60.0;
            }
        }
        hist->Reset();
    }

    gainLogCsv_ << currentTime_ns; rateLogCsv_ << currentTime_ns;
    for (int i = 0; i < 128; ++i) {
        gainLogCsv_ << "," << currentStepPeaks[i];
        rateLogCsv_ << "," << currentStepRates[i];
    }
    gainLogCsv_ << "\n"; rateLogCsv_ << "\n";
    return true;
}

// ヘルパー関数群 (省略なし)
int DetectorAnalyzer::findRightMostPeak(TH1F* hist) {
    int bMin = hist->GetXaxis()->FindBin(PEAK_SEARCH_MIN_TOT);
    for (int i = hist->GetNbinsX(); i >= bMin; --i) if (hist->GetBinContent(i) > 0) return i;
    return 0;
}

void DetectorAnalyzer::determineIntegrationRange(TH1F* hist, int peakBin, int& rMin, int& rMax) {
    if (!hist || peakBin <= 0) { rMin = 0; rMax = 0; return; }
    const int BINS_2US = 100;
    rMax = std::min(hist->GetNbinsX(), peakBin + BINS_2US);
    int foundGapLeftEdge = 0;
    for (int i = peakBin - 1; i >= BINS_2US; --i) {
        bool isGap = true;
        for (int k = 0; k < BINS_2US; ++k) {
            if (hist->GetBinContent(i - k) > 0) { isGap = false; break; }
        }
        if (isGap) { foundGapLeftEdge = i - (BINS_2US - 1); break; }
    }
    if (foundGapLeftEdge == 0) { rMin = 0; rMax = 0; return; }
    rMin = foundGapLeftEdge;
    double widthNs = (double)(rMax - rMin + 1) * 20.0;
    if (widthNs <= 4020.0) { rMin = 0; rMax = 0; }
}

int DetectorAnalyzer::findBestShift(TH1F* hTarget, TH1F* hTemplate, int binMin, int binMax, const std::string& label) {
    double minRes = 1.79769e+308; int bestS = SHIFT_CALC_ERROR; bool found = false;
    const int SEARCH_RANGE = 25; 
    double intTarget = hTarget->Integral(binMin, binMax);
    double intTemplate = hTemplate->Integral(binMin, binMax);
    if (intTarget <= 0 || intTemplate <= 0) return SHIFT_CALC_ERROR;
    double scale = intTemplate / intTarget;

    for (int s = -SEARCH_RANGE; s <= SEARCH_RANGE; ++s) { 
        double res = calculateResidual(hTarget, hTemplate, s, binMin, binMax, scale);
        if (res >= 0 && res < minRes) { minRes = res; bestS = s; found = true; }
    }
    if (!found || std::abs(bestS) == SEARCH_RANGE) return SHIFT_CALC_ERROR;
    return bestS;
}

double DetectorAnalyzer::calculateResidual(TH1F* hTarget, TH1F* hTemplate, int shiftBins, int binMin, int binMax, double scale) {
    int tBinMin = binMin + shiftBins; int tBinMax = binMax + shiftBins;
    if (tBinMin < 1 || tBinMax > hTarget->GetNbinsX()) return -1.0;
    double residualSum = 0.0;
    for (int i = binMin; i <= binMax; ++i) {
        int targetBin = i + shiftBins; 
        double diff = (hTarget->GetBinContent(targetBin) * scale) - hTemplate->GetBinContent(i);
        residualSum += (diff * diff);
    }
    return residualSum;
}

void DetectorAnalyzer::generatePDFReport() {
    if (hRawToTMap.empty()) {
        printLog("[WARN] No data available for PDF report.");
        return;
    }

    printLog("Generating PDF validation report (Summary + Individual Strips)...");
    std::string pdfName = "AnalysisReport.pdf";
    TCanvas* c1 = new TCanvas("c1", "Validation Report", 1200, 800);
    c1->Print((pdfName + "(").c_str()); // PDF作成開始

    // ------------------------------------------------------------------
    // 1. サマリーページ: 検出器タイプ別 ゲイン・ピーク推移 (4枚)
    // ------------------------------------------------------------------
    std::string detNames[] = {"X1", "Y1", "X2", "Y2"};
    c1->SetLogy(0); // 推移グラフは線形スケール
    c1->SetGrid();

    for (int t = 0; t < 4; ++t) {
        c1->Clear();
        TMultiGraph* mg = new TMultiGraph();
        mg->SetTitle(Form("Gain Peak History Summary: %s;Time [min];Peak Position [ToT]", detNames[t].c_str()));

        int graphCount = 0;
        for (int mod = 0; mod < MAX_MODULES; ++mod) {
            for (int ch = 0; ch < MAX_SYS_CH; ++ch) {
                if (!ch_LUT[mod][ch].isValid || ch_LUT[mod][ch].detTypeID != t) continue;
                
                std::pair<int, int> key = {mod, ch};
                if (gainEvolutionGraphs_.count(key) && gainEvolutionGraphs_[key]->GetN() > 0) {
                    TGraph* g = (TGraph*)gainEvolutionGraphs_[key]->Clone();
                    // ストリップ番号に応じて色を変えると見やすい (1-32)
                    int colorIdx = ch_LUT[mod][ch].strip;
                    g->SetLineColor(colorIdx % 9 + 1); 
                    mg->Add(g);
                    graphCount++;
                }
            }
        }

        if (graphCount > 0) {
            mg->Draw("AL"); // A: Axis, L: Line
            c1->Print(pdfName.c_str());
        }
        delete mg;
    }

    // ------------------------------------------------------------------
    // 2. 基本分布 (時間差・イベント幅)
    // ------------------------------------------------------------------
    c1->SetLogy(1);
    if (hDeltaT_Neighbor) { hDeltaT_Neighbor->Draw(); c1->Print(pdfName.c_str()); }
    if (hDeltaT_n8)      { hDeltaT_n8->Draw();      c1->Print(pdfName.c_str()); }
    if (hEventWidth_All)   { hEventWidth_All->Draw();   c1->Print(pdfName.c_str()); }

    // ------------------------------------------------------------------
    // 3. 個別ページ: 対数スペクトル + ゲイン追従範囲(TBox)
    // ------------------------------------------------------------------
    // ストリップ順にソート
    struct SortEntry { int type; int strip; int mod; int ch; };
    std::vector<SortEntry> sortedList;
    for (auto const& [key, rawHist] : hRawToTMap) {
        if (!rawHist) continue;
        auto& cfg = detConfigMap_[key];
        sortedList.push_back({cfg.detTypeID, cfg.strip, key.first, key.second});
    }
    std::sort(sortedList.begin(), sortedList.end(), [](const SortEntry& a, const SortEntry& b) {
        if (a.type != b.type) return a.type < b.type;
        return a.strip < b.strip;
    });

    for (const auto& entry : sortedList) {
        std::pair<int, int> key = {entry.mod, entry.ch};
        TH1F* h = hRawToTMap[key]; 
        
        c1->Clear();
        h->SetTitle(Form("Raw ToT (Log): %s Strip%d (Mod%d Ch%d)", detNames[entry.type].c_str(), entry.strip, entry.mod, entry.ch));
        h->GetXaxis()->SetRangeUser(0, MONITOR_HIST_MAX_TOT);
        h->Draw();

        // ゲイン追従の積分範囲をTBox（赤枠）で重ね描き
        int rMin = rangeMinMap_[key];
        int rMax = rangeMaxMap_[key];
        if (rMin > 0 && rMax > rMin) {
            double xMin = (double)rMin * BIN_WIDTH_NS;
            double xMax = (double)rMax * BIN_WIDTH_NS;
            double yMax = h->GetMaximum() * 1.5;
            TBox* box = new TBox(xMin, 0.5, xMax, yMax);
            box->SetFillStyle(0);
            box->SetLineColor(kRed);
            box->SetLineWidth(2);
            box->Draw("same");
        }
        c1->Print(pdfName.c_str());
    }

    c1->Print((pdfName + ")").c_str()); // PDF作成終了
    delete c1;
    printLog("PDF report generated successfully.");
}


// 他の補助関数も全て維持
void DetectorAnalyzer::printFileTail(const std::string& filePath, long long currentOffset) {
    std::ifstream file(filePath, std::ios::binary);
    if (!file.is_open()) return;
    long long startPos = currentOffset > 64 ? currentOffset - 64 : 0;
    file.seekg(startPos);
    std::vector<unsigned char> buffer(128); 
    file.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
    std::streamsize bytesRead = file.gcount();
    std::cout << "\n========================================================" << std::endl;
    std::cout << " [HEX DUMP] File: " << filePath << std::endl;
    std::cout << " [HEX DUMP] Offset: " << startPos << " - " << (startPos + bytesRead) << std::endl;
    std::cout << "========================================================" << std::endl;
    for (int i = 0; i < bytesRead; ++i) {
        if (i % 16 == 0) std::cout << "\n " << std::setw(8) << std::setfill('0') << std::hex << (startPos + i) << ": ";
        std::cout << std::setw(2) << std::setfill('0') << std::hex << (int)buffer[i] << " ";
        if ((startPos + i) == currentOffset - 1) std::cout << "< "; else std::cout << " ";
    }
    std::cout << "\n========================================================\n" << std::dec << std::endl;
}

void DetectorAnalyzer::calculateEffectiveTime() {
    double liveTimeSec = (double)totalEffectiveTimeNs_ / 1.0e9;
    std::cout << "\n========================================================" << std::endl;
    std::cout << " [RESULT] Total System-wide Live Time: " << std::fixed << std::setprecision(2) << liveTimeSec << " sec" << std::endl;
    double totalElapsed = (double)(currentEventTime_ns_ - globalStartTimestamp_) / 1.0e9;
    std::cout << "          Total Data Span: " << std::fixed << std::setprecision(2) << totalElapsed << " sec" << std::endl;
    if(totalElapsed > 0) {
        std::cout << "          Overall DAQ Efficiency: " << (liveTimeSec / totalElapsed) * 100.0 << " %" << std::endl;
    }
    std::cout << "========================================================\n" << std::endl;
}

double DetectorAnalyzer::calculatePeakAreaCovell(TH1F* hist, double peakPos) {
    if (!hist || peakPos <= 0) return 0.0;
    int binCenter = hist->FindBin(peakPos);
    int halfWidth = COVELL_WINDOW_NS / 20; 
    int binMin = binCenter - halfWidth;
    int binMax = binCenter + halfWidth;
    if (binMin < 1) binMin = 1;
    if (binMax > hist->GetNbinsX()) binMax = hist->GetNbinsX();
    if (binMax <= binMin) return 0.0;
    double grossArea = hist->Integral(binMin, binMax);
    double countLeft  = hist->GetBinContent(binMin);
    double countRight = hist->GetBinContent(binMax);
    double numBins = (double)(binMax - binMin + 1);
    double bgArea = (countLeft + countRight) * numBins / 2.0;
    double netArea = grossArea - bgArea;
    return (netArea < 0) ? 0.0 : netArea;
}

unsigned long long DetectorAnalyzer::getSafeTime(const std::map<int, unsigned long long>& lastTimes) {
    unsigned long long minT = std::numeric_limits<unsigned long long>::max();
    if (lastTimes.empty()) return 0;
    for (auto const& [id, t] : lastTimes) {
        if (t < minT) minT = t;
    }
    return minT;
}

std::string DetectorAnalyzer::getRunSignature(const std::string& fileName) {
    namespace fs = std::filesystem;
    fs::path p(fileName);
    std::string parentDir = p.parent_path().string();
    std::string filenameStr = p.filename().string();
    static std::regex run_re(R"((PSD\d+_\d+))");
    std::smatch match;
    std::string fileRunID = "UnknownID";
    if (std::regex_search(filenameStr, match, run_re)) fileRunID = match[1].str();
    return parentDir + "::" + fileRunID;
}

short DetectorAnalyzer::getRunID(const std::string& prefix) {
    if (runPrefixToId_.count(prefix)) return runPrefixToId_[prefix];
    short newId = (short)runIdToPrefix_.size();
    runIdToPrefix_.push_back(prefix); 
    runPrefixToId_[prefix] = newId;
    return newId;
}

std::string DetectorAnalyzer::parseRunPrefix(const std::string& fullPath) {
    size_t lastSlash = fullPath.find_last_of("/\\");
    std::string fileName = (lastSlash == std::string::npos) ? fullPath : fullPath.substr(lastSlash + 1);
    size_t extPos = fileName.find_last_of(".");
    if (extPos != std::string::npos) fileName = fileName.substr(0, extPos);
    return fileName;
}

void DetectorAnalyzer::loadAndSortFiles(const std::string& directory, std::map<int, std::deque<std::string>>& fileQueues) {
    struct FileInfo { long long date; int runID; int modID; int fileSeq; std::string path;
        bool operator<(const FileInfo& other) const {
            if (date != other.date) return date < other.date;
            if (runID != other.runID) return runID < other.runID;
            return fileSeq < other.fileSeq;
        }
    };
    std::vector<FileInfo> allFiles;
    totalDataSize_ = 0;
    if (!fs::exists(directory)) return;
    printLog("Scanning: Priority=Date -> Run -> Sequence");
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
        if (!entry.is_regular_file() || entry.path().extension() != ".edb") continue;
        std::string fullPath = entry.path().string();
        std::string fileName = entry.path().stem().string(); 
        std::string parentDir = entry.path().parent_path().filename().string(); 
        long long fDate = 0; int fRun = 0;
        size_t pSep = parentDir.find('_');
        if (pSep != std::string::npos) {
            try {
                if (parentDir.substr(0, 3) == "PSD") fRun = std::stoi(parentDir.substr(3, pSep - 3));
                fDate = std::stoll(parentDir.substr(pSep + 1));
            } catch (...) {}
        }
        std::vector<std::string> tokens;
        std::stringstream ss(fileName); std::string item;
        while (std::getline(ss, item, '_')) tokens.push_back(item);
        if (tokens.size() >= 4) {
            try {
                int modID = std::stoi(tokens[tokens.size() - 2]);
                int fileSeq = std::stoi(tokens[tokens.size() - 1]);
                allFiles.push_back({fDate, fRun, modID, fileSeq, fullPath});
                totalDataSize_ += fs::file_size(entry.path());
            } catch (...) {}
        }
    }
    std::sort(allFiles.begin(), allFiles.end());
    for (const auto& info : allFiles) fileQueues[info.modID].push_back(info.path);
    printLog("================================================================");
    printLog("[INFO] Final Sorted File Queue Listing:");
    for (auto const& [mod, q] : fileQueues) {
        printLog("Module " + std::to_string(mod) + ": " + std::to_string(q.size()) + " files found.");
        for (const auto& path : q) printLog("  [Queue Mod " + std::to_string(mod) + "] " + path);
    }
    printLog("================================================================");
}

void DetectorAnalyzer::loadConversionFactorsFromCSV(const std::string& csvFileName, const std::string& detName, int moduleID) {
    std::ifstream file(csvFileName);
    if (!file.is_open()) return;
    int typeID = -1;
    if (detName == "X1") typeID = ID_X1; else if (detName == "Y1") typeID = ID_Y1;
    else if (detName == "X2") typeID = ID_X2; else if (detName == "Y2") typeID = ID_Y2;
    std::string line; std::getline(file, line); 
    if (nameStripToMapKey_.find(detName) == nameStripToMapKey_.end()) nameStripToMapKey_[detName]; 
    activeModuleIDs_.insert(moduleID); 
    struct CSVItem { int localCh; std::pair<int, int> key; std::string label; };
    std::vector<CSVItem> sortedItems;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line); std::string cell; std::vector<std::string> row;
        while (std::getline(ss, cell, ',')) row.push_back(cell);
        if (row.size() < 2) continue;
        try {
            if (row[0].empty() || row[1].empty()) continue;
            int sysCh = std::stoi(row[0]);
            if (detName.find('Y') != std::string::npos) sysCh += Y_CHANNEL_OFFSET; 
            if (sysCh >= MAX_SYS_CH) continue;
            if (row[1].find('N') != std::string::npos || row[1].find('n') != std::string::npos) continue; 
            int localCh = std::stoi(row[1]);
            int pIndex = -1;
            if (row.size() >= 3 && !row[2].empty()) { try { pIndex = std::stoi(row[2]); } catch (...) { pIndex = -1; } }
            std::pair<int, int> keyPair = {moduleID, sysCh};
            detConfigMap_[keyPair] = { typeID, localCh, detName, pIndex };
            if (localCh < 1000) {
                nameStripToMapKey_[detName][localCh] = keyPair;
                std::string label = detName + "_Strip" + std::to_string(localCh);
                sortedItems.push_back({localCh, keyPair, label});
            }
            if (offlineStrips_.count({typeID, localCh})) ignoredMonitorChannels_.insert(keyPair);
            if (moduleID < MAX_MODULES) {
                ch_LUT[moduleID][sysCh].isValid = true;
                ch_LUT[moduleID][sysCh].detTypeID = typeID;
                ch_LUT[moduleID][sysCh].strip = localCh;
                ch_LUT[moduleID][sysCh].pIndex = pIndex;
            }
        } catch (...) { }
    }
    std::sort(sortedItems.begin(), sortedItems.end(), [](const CSVItem& a, const CSVItem& b) { return a.localCh < b.localCh; });
    for (const auto& item : sortedItems) {
        std::string hName = "RawToT_" + item.label;
        std::string hTitle = "Raw ToT: " + item.label;
        if (hRawToTMap.find(item.key) == hRawToTMap.end()) {
            TH1F* h = new TH1F(hName.c_str(), hTitle.c_str(), 100000, 0, 100000);
            h->SetDirectory(nullptr); hRawToTMap[item.key] = h;
        }
        if (gainCheckHists_.find(item.key) == gainCheckHists_.end()) {
            std::string gName = Form("GCheck_%s", item.label.c_str());
            TH1F* gh = new TH1F(gName.c_str(), Form("Gain Check %s", item.label.c_str()), MONITOR_HIST_BINS, 0, MONITOR_HIST_MAX_TOT); 
            gh->SetDirectory(nullptr); gainCheckHists_[item.key] = gh;
        }
        if (gainEvolutionGraphs_.find(item.key) == gainEvolutionGraphs_.end()) {
            TGraph* g = new TGraph();
            g->SetName(Form("GHistory_%s", item.label.c_str()));
            g->SetTitle(Form("Gain Peak History %s;Time [ns];Shift [ns]", item.label.c_str()));
            gainEvolutionGraphs_[item.key] = g;
        }
        if (rateEvolutionGraphs_.find(item.key) == rateEvolutionGraphs_.end()) {
            TGraph* gr = new TGraph();
            gr->SetName(Form("RateHistory_%s", item.label.c_str()));
            gr->SetTitle(Form("Event Rate History %s;Time [ns];Rate [Hz]", item.label.c_str()));
            rateEvolutionGraphs_[item.key] = gr;
        }
        if (item.key.first < MAX_MODULES && item.key.second < MAX_SYS_CH) {
            hRawToT_LUT[item.key.first][item.key.second] = hRawToTMap[item.key];
            hGainCheck_LUT[item.key.first][item.key.second] = gainCheckHists_[item.key];
        }
    }
}

void DetectorAnalyzer::loadOfflineStripList(const std::string& csvFileName) {
    std::ifstream file(csvFileName);
    if (!file.is_open()) return;
    std::string line; std::getline(file, line); 
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line); std::string cell; std::vector<std::string> row;
        while (std::getline(ss, cell, ',')) row.push_back(cell);
        if (row.size() < 3) continue;
        try {
            int detType = std::stoi(row[0]); int stripNo = std::stoi(row[1]); int flag = std::stoi(row[2]);
            if (flag == 1) offlineStrips_.insert({detType, stripNo});
        } catch (...) {}
    }
}

void DetectorAnalyzer::setTimeLimitMinutes(double minutes) {
    if (minutes <= 0) {
        analysisDuration_ns_ = std::numeric_limits<unsigned long long>::max();
        printLog(">>> Time limit set to: UNLIMITED");
    } else {
        analysisDuration_ns_ = static_cast<unsigned long long>(minutes * 60.0 * 1.0e9);
        std::stringstream ss; ss << ">>> Time limit set to: " << std::fixed << std::setprecision(2) << minutes << " minutes";
        printLog(ss.str());
    }
}

void DetectorAnalyzer::setupTree() {
    if (!outputFile_ || !outputFile_->IsOpen()) return;
    outputFile_->cd();
    tree_ = new TTree("tree", "Detector Events");
    tree_->Branch("type", &b_type, "type/B");
    tree_->Branch("strip", &b_strip, "strip/B");
    tree_->Branch("tot", &b_tot, "tot/I");
    tree_->Branch("time", &b_time, "time/L");

    auto initH = [](TH1F*& h, const char* name, const char* title, int bins, double xmin, double xmax) {
        if (h) delete h; h = new TH1F(name, title, bins, xmin, xmax); h->SetDirectory(nullptr);
    };
    initH(hDeltaT_Neighbor, "hDeltaT_Neighbor", "Delta T to Neighbor;Time [ns];Counts", 10000, 0, 10000);
    initH(hDeltaT_n8, "hDeltaT_n8", "Global Delta T (n=8);Time [ns];Counts", 100000, 0, 100000);
    initH(hGlobalPIndexMap, "hGlobalPIndexMap", "Global P-Index Map;P-Index;Counts", 2048, 0, 2048);
    initH(hEventWidth_All, "hEventWidth_All", "Event Width (All);Width [ns];Counts", 1000, 0, 1000);
    initH(hEventWidth_CondA, "hEventWidth_CondA", "Event Width (Cond A);Width [ns];Counts", 1000, 0, 1000);
    initH(hEventWidth_CondB, "hEventWidth_CondB", "Event Width (Cond B);Width [ns];Counts", 1000, 0, 1000);

    if (enableMuonAnalysis_) {
        TDirectory* currentDir = gDirectory;
        idealFile_ = new TFile("ideal_events.root", "RECREATE");
        idealTree_ = new TTree("IdealEvents", "Muon Events");
        b_i_type  = new std::vector<int>(); b_i_strip = new std::vector<int>();
        b_i_tot   = new std::vector<int>(); b_i_time  = new std::vector<double>();
        idealTree_->Branch("runID",   &b_i_runID,   "runID/I");
        idealTree_->Branch("eventID", &b_i_eventID, "eventID/I");
        idealTree_->Branch("nHits",   &b_i_nHits,   "nHits/I");
        idealTree_->Branch("width",   &b_i_width,   "width/D");
        idealTree_->Branch("type",  &b_i_type); idealTree_->Branch("strip", &b_i_strip);
        idealTree_->Branch("tot",   &b_i_tot);  idealTree_->Branch("time",  &b_i_time);
        currentDir->cd();
    }

    for (int mod = 0; mod < MAX_MODULES; ++mod) {
        for (int ch = 0; ch < MAX_SYS_CH; ++ch) {
            if (ch_LUT[mod][ch].isValid) {
                std::pair<int, int> key = {mod, ch};
                initH(hGainCheck_LUT[mod][ch], Form("hGainCheck_M%d_C%d", mod, ch), "Gain Check;ToT;Counts", MONITOR_HIST_BINS, 0, MONITOR_HIST_MAX_TOT);
                gainCheckHists_[key] = hGainCheck_LUT[mod][ch];
                if (gainEvolutionGraphs_.find(key) == gainEvolutionGraphs_.end()) {
                    gainEvolutionGraphs_[key] = new TGraph(); gainEvolutionGraphs_[key]->SetName(Form("gGainEvo_M%d_C%d", mod, ch));
                }
                if (rateEvolutionGraphs_.find(key) == rateEvolutionGraphs_.end()) {
                    rateEvolutionGraphs_[key] = new TGraph(); rateEvolutionGraphs_[key]->SetName(Form("gRateEvo_M%d_C%d", mod, ch));
                }
            }
        }
    }
}

bool DetectorAnalyzer::syncDataStream(std::ifstream& ifs, long long& skippedBytes) {
    skippedBytes = 0; const size_t SCAN_BUFFER_SIZE = 1024 * 1024; 
    std::vector<char> buffer(SCAN_BUFFER_SIZE);
    long long globalOffset = ifs.tellg(); long long initialPos = globalOffset;
    while (ifs) {
        ifs.read(buffer.data(), SCAN_BUFFER_SIZE); size_t bytesRead = ifs.gcount();
        if (bytesRead == 0) break;
        for (size_t i = 0; i < bytesRead; ++i) {
            unsigned char h = static_cast<unsigned char>(buffer[i]);
            if (h == 0x69 || h == 0x6a) {
                bool syncOK = true;
                if (i + 80 <= bytesRead) {
                    for (int k = 1; k < 10; ++k) {
                        unsigned char nextH = static_cast<unsigned char>(buffer[i + (k * 8)]);
                        if (nextH != 0x69 && nextH != 0x6a) { syncOK = false; break; }
                    }
                } else {
                    long long candidatePos = globalOffset + i; long long nextChunkPos = globalOffset + bytesRead; 
                    for (int k = 1; k < 10; ++k) {
                         ifs.seekg(candidatePos + (k * 8), std::ios::beg); char nextH;
                         if (!ifs.read(&nextH, 1)) { syncOK = false; break; }
                         if ((unsigned char)nextH != 0x69 && (unsigned char)nextH != 0x6a) { syncOK = false; break; }
                    }
                    if (!syncOK) { ifs.clear(); ifs.seekg(nextChunkPos, std::ios::beg); }
                }
                if (syncOK) {
                    long long finalPos = globalOffset + i; ifs.seekg(finalPos, std::ios::beg);
                    skippedBytes = finalPos - initialPos; return true;
                }
            }
        }
        globalOffset += bytesRead;
    }
    skippedBytes = globalOffset - initialPos; return false;
}