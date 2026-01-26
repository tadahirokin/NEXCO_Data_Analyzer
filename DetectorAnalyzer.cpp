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

    for (int m : activeModuleIDs_) {
    lastT0_[m] = 0;
    hasDataAligned_[m] = false;
    hasFoundT0_[m] = false;
    
    // ★ これを追加
    lastRawT0_[m] = 0;
    moduleOffsets_[m] = 0;
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


static std::map<int, bool> g_allowNextT0Jump; 






void DetectorAnalyzer::printSearchStatus() {
    double liveTimeMin = (double)totalEffectiveTimeNs_ / 1.0e9 / 60.0;
    double progress = (analysisDuration_ns_ > 0 && analysisDuration_ns_ != ULLONG_MAX) 
                      ? (double)totalEffectiveTimeNs_ / (double)analysisDuration_ns_ * 100.0 : 0.0;
    std::cout << "\r\033[K" << "Progress: [" << std::fixed << std::setprecision(1) << progress << "%] "
              << "Live Time: " << liveTimeMin << " min" << std::flush;
}



void DetectorAnalyzer::printSurvivalReport() {
    std::cout << "\n" << std::string(70, '=') << "\n [SURVIVAL REPORT at 5.0 min]\n" << std::string(70, '-') << std::endl;
    std::string detNames[] = {"X1", "Y1", "X2", "Y2"};
    for (int t = 0; t < 4; ++t) {
        int alive = 0, total = 0;
        std::vector<int> silent;
        for (auto const& [key, h] : hRawToTMap) {
            if (ch_LUT[key.first][key.second].detTypeID == t) {
                total++;
                if (aliveStatus_[key]) alive++; else silent.push_back(ch_LUT[key.first][key.second].strip);
            }
        }
        std::cout << " " << detNames[t] << ": " << alive << " / " << total << " Alive. ";
        if (!silent.empty()) {
            std::cout << "Silent Strips: ";
            std::sort(silent.begin(), silent.end());
            for(size_t j=0; j<silent.size(); ++j) std::cout << silent[j] << (j==silent.size()-1?"":",");
        }
        std::cout << std::endl;
    }
    std::cout << std::string(70, '=') << "\n" << std::endl;
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


// ----------------------------------------------------------------------------
// 1. readEventsFromFile (シーク方式による復元版)
//    現在の最新コードの同名関数を、このコードで上書きしてください。
//    これにより、バッファ境界バグを回避しつつ、ヒストグラム機能等は維持されます。
// ----------------------------------------------------------------------------
std::pair<unsigned long long, bool> DetectorAnalyzer::readEventsFromFile(const std::string& fileName, int modID, std::vector<Event>& rawEvents, long long& offset) {
    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs.is_open()) return {lastT0_[modID], true};

    // ファイルサイズと終端の確認
    ifs.seekg(0, std::ios::end);
    long long fileSize = static_cast<long long>(ifs.tellg());
    if (offset >= fileSize) return {lastT0_[modID], true};

    // 64MB一括読み込み
    const size_t CHUNK_SIZE = 64 * 1024 * 1024;
    ifs.seekg(offset, std::ios::beg);
    std::vector<unsigned char> buf(CHUNK_SIZE);
    ifs.read(reinterpret_cast<char*>(buf.data()), CHUNK_SIZE);
    size_t count = static_cast<size_t>(ifs.gcount());
    bool isEOF = (offset + static_cast<long long>(count) >= fileSize);

    if (count < 80 && !hasDataAligned_[modID]) {
        offset += static_cast<long long>(count);
        return {lastT0_[modID], isEOF};
    }

    // --- Phase 1: 物理アライメント探索 (10パケット連続ヘッダ確認) ---
    if (!hasDataAligned_[modID]) {
        for (size_t i = 0; i + 80 <= count; ++i) {
            if (buf[i] == 0x69 || buf[i] == 0x6a) {
                bool match = true;
                for (int k = 1; k < 10; ++k) {
                    if (buf[i + k * 8] != 0x69 && buf[i + k * 8] != 0x6a) {
                        match = false; break;
                    }
                }
                if (match) {
                    hasDataAligned_[modID] = true;
                    offset += static_cast<long long>(i);
                    // 同期直後はオフセットもリセット（新しいファイルの基準）
                    moduleOffsets_[modID] = 0; 
                    return {0, false}; 
                }
            }
        }
        size_t advance = (count > 80) ? (count - 80) : count;
        offset += static_cast<long long>(advance);
        return {lastT0_[modID], isEOF};
    }

    // --- Phase 2: 高速デコード ---
    size_t i = 0;
    unsigned long long ultimateT0_corrected = lastT0_[modID];    
    unsigned long long currentBaseTime_corrected = lastT0_[modID];

    while (i + 8 <= count) {
        unsigned char header = buf[i];

        if (header != 0x69 && header != 0x6a) {
            hasDataAligned_[modID] = false; 
            hasFoundT0_[modID] = false; 
            offset += (static_cast<long long>(i) + 1); 
            return {ultimateT0_corrected, isEOF};
        }

        // [T0 パケット処理]
        if (header == 0x69) { 
            unsigned char* p = &buf[i + 1];
            // --- 生時刻のデコード (raw_t) ---
            unsigned long long raw_t = ((unsigned long long)p[0] << 24 | (unsigned long long)p[1] << 16 | 
                                       (unsigned long long)p[2] << 8 | (unsigned long long)p[3]) * 1000000000ULL + 
                                       ((unsigned long long)p[4] << 2 | (unsigned long long)(p[5] & 0xC0) >> 6) * 1000000ULL + 
                                       ((unsigned long long)(p[5] & 0x3F) << 2 | (unsigned long long)(p[6] & 0xC0) >> 6);

            // --- ワープ検知 & 10パケット検証ロジック (シーク方式) ---
            if (hasFoundT0_[modID]) {
                // lastRawT0_[modID] との差分を確認 (1msリズムから外れたか)
                long long diff = (long long)raw_t - (long long)lastRawT0_[modID];
                
                if (std::abs(diff - 1000000LL) > 500000LL) {
                    // ワープの疑いあり。未来を覗き見して検証 (ディスクシークを使用)
                    std::streampos currentFilePos = ifs.tellg(); // チャンク読み込み後の位置
                    
                    // 実際の位置 = offset + 現在のバッファ内位置(i) + 次のパケットへ(8)
                    ifs.seekg(offset + i + 8, std::ios::beg); 
                    
                    int validT0Count = 0;
                    unsigned long long tempPrevRaw = raw_t;
                    unsigned char v[8];
                    
                    // 次のT0を10個見つけてリズムを確認
                    while (validT0Count < 10 && ifs.read((char*)v, 8)) {
                        if (v[0] == 0x69) {
                            unsigned long long vt = ((unsigned long long)v[1] << 24 | (unsigned long long)v[2] << 16 | 
                                                    (unsigned long long)v[3] << 8 | (unsigned long long)v[4]) * 1000000000ULL + 
                                                    ((unsigned long long)v[5] << 2 | (unsigned long long)(v[6] & 0xC0) >> 6) * 1000000ULL + 
                                                    ((unsigned long long)(v[6] & 0x3F) << 2 | (unsigned long long)(v[7] & 0xC0) >> 6);
                            
                            if (std::abs((long long)(vt - tempPrevRaw) - 1000000LL) < 100000LL) {
                                validT0Count++;
                                tempPrevRaw = vt;
                            } else { break; }
                        }
                    }
                    // 元の位置（チャンク末端）に戻す (これを忘れるとファイル読み込みが狂うが、ここは正しい)
                    ifs.seekg(currentFilePos, std::ios::beg); 
                    ifs.clear();

                    if (validT0Count >= 10) {
                        // 本物のワープ確定。新しいオフセットを算出
                        long long jump = (long long)raw_t - ((long long)lastRawT0_[modID] + 1000000LL);
                        moduleOffsets_[modID] += jump;
                        
                        // ログ出力 (デバッグ用)
                        // std::cout << "[WARP] Mod" << modID << " Jump corrected. Offset=" << moduleOffsets_[modID] << std::endl;
                    }
                }
            }

            // --- 補正後時刻の計算 ---
            unsigned long long corrected_t = raw_t - moduleOffsets_[modID];
            
            i += 8;
            if (corrected_t <= lastT0_[modID] && lastT0_[modID] != 0) continue; // 重複排除

            hasFoundT0_[modID] = true;
            lastRawT0_[modID] = raw_t;       // 次回判定用に「生」を保存
            lastT0_[modID] = corrected_t;    // クラス全体には「補正後」を報告
            ultimateT0_corrected = corrected_t;
            currentBaseTime_corrected = corrected_t;
            continue;
        }

        // [Event パケット処理]
        if (header == 0x6a) { 
            if (hasFoundT0_[modID] && currentBaseTime_corrected > 0) {
                unsigned char* p = &buf[i + 1];
                unsigned int tof = ((unsigned int)p[0] << 16 | (unsigned int)p[1] << 8 | (unsigned int)p[2]);
                unsigned int pw  = ((unsigned int)p[3] << 12 | (unsigned int)p[4] << 4 | (unsigned int)(p[5] & 0xF0) >> 4);
                int pixel = (int)p[6];
                
                // すでに補正済みの currentBaseTime_corrected を使うので、evTime も自動的に補正される
                unsigned long long evTime = currentBaseTime_corrected + static_cast<unsigned long long>(tof) - static_cast<unsigned long long>(pw);
                rawEvents.push_back({0, pixel, modID, evTime, (int)pw, pixel});
            }
            i += 8;
        }
    }

    offset += static_cast<long long>(count);
    return {ultimateT0_corrected, isEOF};
}

// ----------------------------------------------------------------------------
// 3. determineIntegrationRange (修正版)
//    左側の空白区間(2μs)を積分範囲に「含める」ように論理修正
// ----------------------------------------------------------------------------
void DetectorAnalyzer::determineIntegrationRange(TH1F* hist, int peakBin, int& outMinBin, int& outMaxBin) {
    int nBins = hist->GetNbinsX();
    const int SILENCE_BINS = 100; // 2µs @ 20ns/bin
    
    // --- 左側の探索 ---
    outMinBin = 1; 
    int zeroCount = 0;
    
    // ピークから左に向かって走査
    for (int i = peakBin; i > 1; --i) {
        if (hist->GetBinContent(i) <= 0) zeroCount++; 
        else zeroCount = 0;
        
        // 2μs分の空白が見つかったら
        if (zeroCount >= SILENCE_BINS) { 
            // ★修正: i + zeroCount (信号の際) ではなく、
            // i (空白の左端) を採用して、マージンを確保する
            outMinBin = i; 
            break; 
        }
    }
    
    // --- 右側の探索 ---
    outMaxBin = nBins; 
    zeroCount = 0;
    
    // ピークから右に向かって走査
    for (int i = peakBin; i < nBins; ++i) {
        if (hist->GetBinContent(i) <= 0) zeroCount++; 
        else zeroCount = 0;
        
        if (zeroCount >= SILENCE_BINS) { 
            outMaxBin = i; // こちらは空白を含んだ位置でOK
            break; 
        }
    }
}


// ----------------------------------------------------------------------------
// 2. processChunk
//    絶対時刻監視 & 10分トリガー
// ----------------------------------------------------------------------------
bool DetectorAnalyzer::processChunk(const std::vector<Event>& sortedEvents) {
    if (sortedEvents.empty()) return true;

    for (const auto& e : sortedEvents) {
        if (currentRunStartTime_ == 0 && e.tot >= 1000) {
            currentRunStartTime_ = e.eventTime_ns;
            if (lastGlobalAnalysisTime_ns_ == 0) {
                lastGlobalAnalysisTime_ns_ = e.eventTime_ns;
            }
            for(int mid : activeModuleIDs_) moduleAliveTime_[mid] = e.eventTime_ns;
        }
        currentEventTime_ns_ = e.eventTime_ns;
        moduleAliveTime_[e.module] = e.eventTime_ns;

        // 10分 (GAIN_ANALYSIS_WINDOW_NS) 以上経過したら解析実行
        if (e.eventTime_ns >= lastGlobalAnalysisTime_ns_ + GAIN_ANALYSIS_WINDOW_NS) {
            analyzeGainShift(e.eventTime_ns);
            lastGlobalAnalysisTime_ns_ = e.eventTime_ns;
        }

        if (e.module >= 0 && e.module < MAX_MODULES && e.sysCh >= 0 && e.sysCh < MAX_SYS_CH) {
            TH1F* h = hRawToT_LUT[e.module][e.sysCh];
            if (h != nullptr) {
                h->Fill((double)e.tot); 
                if (hGainCheck_LUT[e.module][e.sysCh]) hGainCheck_LUT[e.module][e.sysCh]->Fill((double)e.tot);
                aliveStatus_[{e.module, e.sysCh}] = true;

                if (tree_) {
                    b_type  = static_cast<Char_t>(ch_LUT[e.module][e.sysCh].detTypeID);
                    b_strip = static_cast<Char_t>(ch_LUT[e.module][e.sysCh].strip);
                    b_tot   = static_cast<Int_t>(e.tot);
                    b_time  = static_cast<Long64_t>(e.eventTime_ns); 
                    tree_->Fill();
                }
            }
        }
        
        if (enableMuonAnalysis_) {
            analysisBuffer_.push_back(e);
        }
    }
    return true;
}

// ----------------------------------------------------------------------------
// 4. processBinaryFiles (完全版)
//    ・TotalEffectiveTime廃止 -> ランごとのLiveTime積算方式へ変更
//    ・lastRawT0_, moduleOffsets_ の初期化漏れ修正
//    ・ゲイン解析トリガーを SafeTime (絶対時刻) ベースに変更
// ----------------------------------------------------------------------------
void DetectorAnalyzer::processBinaryFiles(std::map<int, std::deque<std::string>>& fileQueues) {
    const size_t PER_MOD_CAP = 10000000;
    std::map<int, std::deque<Event>> moduleBuffers;
    std::map<int, unsigned long long> moduleLastTime;
    auto lastUIDraw = std::chrono::system_clock::now();

    // --- 新・時間管理変数 ---
    unsigned long long accumulatedLiveTime_ns = 0;     // 完了した過去のランの合計ライブタイム
    unsigned long long currentRunStartSafeTime = 0;    // 現在のランの開始T0 (SafeTime)
    unsigned long long lastLoopSafeTime = 0;           // 直前ループのT0 (ラン切り替わり/逆転検知用)
    bool isFirstSafeTime = true;                       // 解析開始直後のフラグ

    // アクティブなモジュールIDを特定
    activeModuleIDs_.clear();
    for(auto const& [m, q] : fileQueues) if(!q.empty()) activeModuleIDs_.insert(m);
    
    // 変数初期化 (ループ内で確実にリセット)
    for (int m : activeModuleIDs_) {
        lastT0_[m] = 0; 
        hasDataAligned_[m] = false; 
        hasFoundT0_[m] = false;
        moduleLastTime[m] = 0; 
        currentFileOffsets_[m] = 0;
        
        // ★修正: 初期化漏れにより異常値が出ていた変数をここで確実にゼロ化
        lastRawT0_[m] = 0;     
        moduleOffsets_[m] = 0; 
    }

    std::string currentRunSignature = "";
    for(auto const& [m, q] : fileQueues) if(!q.empty()) { currentRunSignature = getRunSignature(q.front()); break; }

    while (true) {
        // --- 1. データ読み込みフェーズ ---
        int lagMod = -1; unsigned long long minT = 0xFFFFFFFFFFFFFFFFULL;
        for (int m : activeModuleIDs_) {
            if (fileQueues[m].empty() && moduleBuffers[m].empty()) continue;
            if (moduleLastTime[m] < minT) { minT = moduleLastTime[m]; lagMod = m; }
        }
        if (lagMod == -1) break; // 全モジュール完了

        // バッファ補充
        if (lagMod != -1 && moduleBuffers[lagMod].size() < PER_MOD_CAP && !fileQueues[lagMod].empty()) {
            if (getRunSignature(fileQueues[lagMod].front()) == currentRunSignature) {
                std::vector<Event> temp;
                auto res = readEventsFromFile(fileQueues[lagMod].front(), lagMod, temp, currentFileOffsets_[lagMod]);
                if (!temp.empty()) {
                    moduleBuffers[lagMod].insert(moduleBuffers[lagMod].end(), std::make_move_iterator(temp.begin()), std::make_move_iterator(temp.end()));
                    processedDataSize_ += (64 * 1024 * 1024);
                }
                moduleLastTime[lagMod] = res.first;
                
                if (res.second) { // EOF
                    printLog("[FILE DONE] " + fileQueues[lagMod].front());
                    fileQueues[lagMod].pop_front(); 
                    currentFileOffsets_[lagMod] = 0;
                    // ファイルが変わっても hasFoundT0_ は維持する (連続性のため)
                    // ただしアライメントは再確認させる
                    hasDataAligned_[lagMod] = false; 
                }
            }
        }

        // --- 2. 同期時刻 (SafeTime) の取得 ---
        unsigned long long safeTime = getSafeTime(moduleLastTime);

        if (safeTime > 0) {
            // --- 3. ライブタイム計算ロジック (ラン積算方式) ---
            
            // 初回検出時、または時間が巻き戻った場合 (新しいランの開始)
            if (isFirstSafeTime || safeTime < lastLoopSafeTime) {
                if (!isFirstSafeTime) {
                    // 直前のランの時間を確定して積算
                    unsigned long long runDuration = (lastLoopSafeTime > currentRunStartSafeTime) 
                                                   ? (lastLoopSafeTime - currentRunStartSafeTime) : 0;
                    accumulatedLiveTime_ns += runDuration;
                    
                    // ゲイン解析用の次回トリガーをリセット (新しいランの開始に合わせて)
                    lastGlobalAnalysisTime_ns_ = safeTime; 
                }
                currentRunStartSafeTime = safeTime;
                isFirstSafeTime = false;
            }

            // 現在進行中のランの経過時間
            unsigned long long currentRunDuration = (safeTime >= currentRunStartSafeTime) 
                                                  ? (safeTime - currentRunStartSafeTime) : 0;
            
            // 全体の有効ライブタイム (表示・制限用)
            finalTotalTimeNs_ = accumulatedLiveTime_ns + currentRunDuration;

            // 直前時刻を更新
            lastLoopSafeTime = safeTime;


            // --- 4. 終了判定 ---
            if (analysisDuration_ns_ != std::numeric_limits<unsigned long long>::max()) {
                if (finalTotalTimeNs_ >= analysisDuration_ns_) break;
            }


            // --- 5. ゲイン解析 (SafeTime 絶対時刻ベース) ---
            if (lastGlobalAnalysisTime_ns_ == 0) lastGlobalAnalysisTime_ns_ = safeTime;

            // 時間が逆転していない、かつ 前回から一定時間経過していたら実行
            if (safeTime >= lastGlobalAnalysisTime_ns_) {
                if (safeTime - lastGlobalAnalysisTime_ns_ >= GAIN_ANALYSIS_WINDOW_NS) {
                    analyzeGainShift(safeTime); // 絶対時刻を渡す
                    lastGlobalAnalysisTime_ns_ = safeTime;
                }
            } else {
                // ラン切り替わり等で SafeTime が若返った場合、基準時刻を現在にリセット
                lastGlobalAnalysisTime_ns_ = safeTime;
            }


            // --- 6. イベント処理 (Chunking) ---
            std::vector<Event> chunk;
            for (int m : activeModuleIDs_) {
                while (!moduleBuffers[m].empty() && moduleBuffers[m].front().eventTime_ns <= safeTime) {
                    chunk.push_back(std::move(moduleBuffers[m].front()));
                    moduleBuffers[m].pop_front();
                }
            }
            if (!chunk.empty()) {
                std::sort(chunk.begin(), chunk.end(), [](const Event& a, const Event& b) { return a.eventTime_ns < b.eventTime_ns; });
                processChunk(chunk);
            }
        }

        // --- 7. UI更新 (200ms) ---
        auto now = std::chrono::system_clock::now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastUIDraw).count() > 200) {
            lastUIDraw = now;
            
            double progress = 0.0;
            if (analysisDuration_ns_ != std::numeric_limits<unsigned long long>::max()) {
                progress = (double)finalTotalTimeNs_ / (double)analysisDuration_ns_ * 100.0;
            } else if (totalDataSize_ > 0) {
                progress = (double)processedDataSize_ / totalDataSize_ * 100.0;
            }
            if (progress > 100.0) progress = 100.0;

            std::stringstream ss;
            ss << "\r\033[K" 
               << "Progress: [" << std::fixed << std::setprecision(1) << progress << "%] "
               << "Live: " << std::fixed << std::setprecision(2) << (double)finalTotalTimeNs_/1e9/60.0 << " min "
               << "(SafeTime: " << safeTime << ")\n";
            
            int lineCount = 1;
            for (int m : activeModuleIDs_) {
                auto getBar = [&](int mID) {
                    int w = 15; 
                    double r = (double)moduleBuffers[mID].size() / PER_MOD_CAP;
                    if (r > 1.0) r = 1.0;
                    int f = (int)(r * w);
                    std::string b = "["; for(int k=0; k<w; ++k) b += (k < f ? "#" : "."); return b + "]";
                };
                auto getFile = [&](int mID) {
                    if (fileQueues[mID].empty()) return std::string("(No File)");
                    return std::filesystem::path(fileQueues[mID].front()).filename().string();
                };
                ss << "\033[K" << " M" << m << ": " << std::left << std::setw(28) << getFile(m).substr(0, 27) << " " << getBar(m) << "\n";
                lineCount++;
            }
            ss << "\033[" << lineCount << "A"; 
            std::cout << ss.str() << std::flush;
        }
    }
    std::cout << "\033[" << (activeModuleIDs_.size() + 1) << "B" << std::endl;
}

// ----------------------------------------------------------------------------
// ゲイン解析 (SafeTime対応版)
// 引数 currentTime_ns は絶対時刻(T0)を受け取る前提
// ----------------------------------------------------------------------------
bool DetectorAnalyzer::analyzeGainShift(unsigned long long currentTime_ns) {
    // グラフのX軸用 (分単位)
    double timeMin = (double)currentTime_ns / 60.0e9;

    // レート計算用の経過時間 (秒)
    // 前回実行時からの差分を取る。SafeTimeベースなので正確。
    double durationSec = 0.0;
    
    // 初回、またはラン切り替わりで時間が戻った場合は、デフォルト値(600s)か 1.0s を使う
    if (lastGlobalAnalysisTime_ns_ > 0 && currentTime_ns > lastGlobalAnalysisTime_ns_) {
        durationSec = (double)(currentTime_ns - lastGlobalAnalysisTime_ns_) / 1.0e9;
    } else {
        durationSec = 600.0; // デフォルト (10分間隔前提)
    }
    if (durationSec < 1.0) durationSec = 1.0; // 除算エラー防止

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

    // 全チャネルスキャン
    for (auto const& [key, hist] : gainCheckHists_) {
        auto& cfg = detConfigMap_[key];
        if (!hist) continue;

        // 1. テンプレート未作成なら作成 (初回)
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
        // 2. テンプレートありならマッチング
        else if (rangeMinMap_[key] > 0) {
            TH1F* tHist = templateHists_[key];
            int deltaBin = findBestShift(hist, tHist, rangeMinMap_[key], rangeMaxMap_[key]);
            
            if (deltaBin != SHIFT_CALC_ERROR) {
                rangeMinMap_[key] += deltaBin;
                rangeMaxMap_[key] += deltaBin;
                cumulativeShiftNsMap_[key] += (deltaBin * 20.0);

                double absPeakPos = initialPeakPosMap_[key] + cumulativeShiftNsMap_[key];
                if (!gainEvolutionGraphs_[key]) gainEvolutionGraphs_[key] = new TGraph();
                gainEvolutionGraphs_[key]->SetPoint(gainEvolutionGraphs_[key]->GetN(), timeMin, absPeakPos);
                
                // テンプレート更新 (追従)
                tHist->Reset();
                tHist->Add(hist);
            }
        }

        // 3. データ記録 (ピーク位置 & CPM)
        if (rangeMinMap_[key] > 0) {
            int idx = cfg.detTypeID * 32 + (cfg.strip - 1);
            if (idx >= 0 && idx < 128) {
                currentStepPeaks[idx] = initialPeakPosMap_[key] + cumulativeShiftNsMap_[key];
                
                // 積分範囲内のカウント数
                double counts = hist->Integral(rangeMinMap_[key], rangeMaxMap_[key]);
                // CPM計算: (Counts / DurationSec) * 60
                currentStepRates[idx] = (counts / durationSec) * 60.0;
            }
        }
        
        // 次の区間のためにヒストグラムをリセット
        hist->Reset();
    }

    // CSV出力
    gainLogCsv_ << currentTime_ns; rateLogCsv_ << currentTime_ns;
    for (int i = 0; i < 128; ++i) {
        gainLogCsv_ << "," << currentStepPeaks[i];
        rateLogCsv_ << "," << currentStepRates[i];
    }
    gainLogCsv_ << "\n"; rateLogCsv_ << "\n";
    return true;
}

// ----------------------------------------------------------------------------
// 結果表示 (LiveTime表示の更新)
// ----------------------------------------------------------------------------
void DetectorAnalyzer::calculateEffectiveTime() {
    double liveTimeSec = (double)finalTotalTimeNs_ / 1.0e9;
    
    // 設定時間
    double targetSec = (analysisDuration_ns_ != std::numeric_limits<unsigned long long>::max()) 
                       ? (double)analysisDuration_ns_ / 1.0e9 : 0.0;

    std::cout << "\n========================================================" << std::endl;
    std::cout << " [RESULT] Total System-wide Live Time: " << std::fixed << std::setprecision(2) << liveTimeSec << " sec" << std::endl;
    std::cout << "          (Accumulated from synchronized T0 intervals)" << std::endl;
    
    if (targetSec > 0) {
        double ratio = (liveTimeSec / targetSec) * 100.0;
        std::cout << "          Target Completion Ratio: " << std::fixed << std::setprecision(2) << ratio << " % (Limit: " << targetSec << "s)" << std::endl;
    }
    std::cout << "========================================================\n" << std::endl;
}


// ヘルパー関数群 (省略なし)
int DetectorAnalyzer::findRightMostPeak(TH1F* hist) {
    int bMin = hist->GetXaxis()->FindBin(PEAK_SEARCH_MIN_TOT);
    for (int i = hist->GetNbinsX(); i >= bMin; --i) if (hist->GetBinContent(i) > 0) return i;
    return 0;
}


int DetectorAnalyzer::findBestShift(TH1F* hTarget, TH1F* hTemplate, int binMin, int binMax) {
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


// ----------------------------------------------------------------------------
// 6. generatePDFReport (修正版)
//    機能追加: 参照光ピークの移動軌跡を、時間経過＝高さ として重ね描きする
// ----------------------------------------------------------------------------
void DetectorAnalyzer::generatePDFReport() {
    if (hRawToTMap.empty()) {
        printLog("[WARN] No data available for PDF report.");
        return;
    }

    printLog("Generating PDF validation report (128 Strips + Summaries)...");
    std::string pdfName = "AnalysisReport.pdf";
    TCanvas* c1 = new TCanvas("c1", "Final Report", 1200, 800);
    c1->Print((pdfName + "(").c_str()); 

    std::string detNames[] = {"X1", "Y1", "X2", "Y2"};

    // 1. 個別ページ (128ページ)
    c1->SetLogy(1); // ログスケール
    
    // ソート用構造体
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
        h->SetTitle(Form("%s Strip%d (Mod%d Ch%d);ToT;Counts (Overlay: Peak Drift)", detNames[entry.type].c_str(), entry.strip, entry.mod, entry.ch));
        h->GetXaxis()->SetRangeUser(0, MONITOR_HIST_MAX_TOT);
        h->SetLineColor(kBlue+1);
        h->Draw();

        double yMax = h->GetMaximum() * 1.5;
        double yMin = 0.5; // Logスケールの下限目安

        // --- A. 積分範囲のボックス (赤枠) ---
        int rMin = rangeMinMap_[key];
        int rMax = rangeMaxMap_[key];
        if (rMin > 0 && rMax > rMin) {
            double xMin = (double)rMin * 20.0; 
            double xMax = (double)rMax * 20.0;
            TBox* box = new TBox(xMin, yMin, xMax, yMax);
            box->SetFillStyle(0);
            box->SetLineColor(kRed);
            box->SetLineWidth(2);
            box->Draw("same");
        }

        // --- B. ピーク移動の軌跡 (緑の点) ---
        // 時間経過をY軸の高さ（対数）に変換してプロット
        if (gainEvolutionGraphs_.count(key)) {
            TGraph* gEvol = gainEvolutionGraphs_[key];
            if (gEvol->GetN() > 0) {
                TGraph* gVis = new TGraph();
                
                // グラフの最後の時間が「最高到達点」になるように正規化
                double tMax = 0;
                double x_dummy, t_last;
                gEvol->GetPoint(gEvol->GetN()-1, t_last, x_dummy);
                tMax = t_last; 
                if (tMax <= 0) tMax = 1.0; // 安全策

                double logMin = std::log10(yMin);
                double logMax = std::log10(yMax);

                for(int i=0; i < gEvol->GetN(); ++i) {
                    double t_ns, pos_tot;
                    gEvol->GetPoint(i, t_ns, pos_tot);
                    
                    // 時間(0 -> tMax) を Y軸(yMin -> yMax) に対数でマッピング
                    double ratio = t_ns / tMax;
                    double logY = logMin + (logMax - logMin) * ratio;
                    double plotY = std::pow(10.0, logY);

                    gVis->SetPoint(i, pos_tot, plotY);
                }
                
                gVis->SetMarkerStyle(20); // 塗りつぶし円
                gVis->SetMarkerSize(0.6);
                gVis->SetMarkerColor(kGreen+2); // 緑色
                gVis->Draw("P same");
            }
        }

        c1->Print(pdfName.c_str());
    }

    // 2. 時間差・空間分布
    c1->SetLogy(1);
    if (hDeltaT_Neighbor) { hDeltaT_Neighbor->GetXaxis()->SetRangeUser(0, 1000000); hDeltaT_Neighbor->Draw(); c1->Print(pdfName.c_str()); }
    if (hDeltaT_n8) { hDeltaT_n8->GetXaxis()->SetRangeUser(0, 1000000); hDeltaT_n8->Draw(); c1->Print(pdfName.c_str()); }
    c1->SetLogy(0);
    if (hGlobalPIndexMap) { hGlobalPIndexMap->Draw(); c1->Print(pdfName.c_str()); }

    // 3. サマリーグラフ
    for (int t = 0; t < 4; ++t) {
        c1->Clear(); c1->SetGrid();
        TMultiGraph* mgG = new TMultiGraph();
        mgG->SetTitle(Form("Gain Peak History: %s;Time [min];Peak [ToT]", detNames[t].c_str()));
        for (auto const& [key, g] : gainEvolutionGraphs_) {
            if (detConfigMap_[key].detTypeID == t && g->GetN() > 0) {
                TGraph* gc = (TGraph*)g->Clone();
                gc->SetLineColor(detConfigMap_[key].strip % 9 + 1);
                mgG->Add(gc);
            }
        }
        if (mgG->GetListOfGraphs()) { mgG->Draw("AL"); c1->Print(pdfName.c_str()); }
    }

    for (int t = 0; t < 4; ++t) {
        c1->Clear(); c1->SetGrid();
        TMultiGraph* mgR = new TMultiGraph();
        mgR->SetTitle(Form("Rate History (CPM): %s;Time [min];Rate [CPM]", detNames[t].c_str()));
        for (auto const& [key, g] : rateEvolutionGraphs_) {
            if (detConfigMap_[key].detTypeID == t && g->GetN() > 0) {
                TGraph* rc = (TGraph*)g->Clone();
                rc->SetLineColor(detConfigMap_[key].strip % 9 + 1);
                mgR->Add(rc);
            }
        }
        if (mgR->GetListOfGraphs()) { mgR->Draw("AL"); c1->Print(pdfName.c_str()); }
    }

    c1->Print((pdfName + ")").c_str()); 
    delete c1;
    printLog("Full PDF report generated successfully.");
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


// --- クラスのメンバ関数としての printLog ---
void DetectorAnalyzer::printLog(const std::string& msg) {
    std::cout << "\r\033[K" << std::flush;
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm* now_tm = std::localtime(&now_c);
    std::cout << "[" << std::put_time(now_tm, "%Y-%m-%d %H:%M:%S") << "] [INFO] " << msg << std::endl;
}

// ----------------------------------------------------------------------------
// 1. setupTree: ヒストグラムの先行生成 (順序保証・高解像度)
// ----------------------------------------------------------------------------
void DetectorAnalyzer::setupTree() {
    if (!outputFile_ || !outputFile_->IsOpen()) return;
    outputFile_->cd();

    // Tree構造の定義 (Legacy Codeの型定義に準拠)
    if (tree_) delete tree_;
    tree_ = new TTree("tree", "Detector Events");
    tree_->Branch("type",  &b_type,  "type/B");   // Char_t
    tree_->Branch("strip", &b_strip, "strip/B");  // Char_t
    tree_->Branch("tot",   &b_tot,   "tot/I");    // Int_t
    tree_->Branch("time",  &b_time,  "time/L");   // Long64_t

    // --- ヒストグラムの先行生成 ---
    // X1(1-32) -> Y1(1-32) -> X2(1-32) -> Y2(1-32) の順で生成し、ROOTファイル内の順序を固定する
    const char* detNames[] = {"X1", "Y1", "X2", "Y2"};
    
    // 全128チャンネル分を確保
    for (int type = 0; type < 4; ++type) {
        for (int strip = 1; strip <= 32; ++strip) {
            // Global Index: 0-127 (配列管理用)
            int gIdx = (type * 32) + (strip - 1);
            
            // 名前とタイトル: 物理名を明記
            // 例: hRaw_X1_Strip01
            std::string name = Form("hRaw_%s_Strip%02d", detNames[type], strip);
            std::string title = Form("%s Strip %02d;ToT [ns];Counts", detNames[type], strip);

            // 既存があれば削除 (二重生成防止)
            if (hRawToT_Global[gIdx]) delete hRawToT_Global[gIdx];
            
            // ★重要: 100,000 bins, 0 - 100,000 ns (1ns/bin)
            hRawToT_Global[gIdx] = new TH1F(name.c_str(), title.c_str(), 100000, 0, 100000);
            
            // デフォルトではディレクトリに関連付けない（書き込み時に手動で制御する場合）
            // ただし今回は最後に一括Writeするので、ディレクトリに関連付けても良いが、
            // あかりさんの流儀に合わせて管理します。
            
            // ゲインチェック用ヒストグラムも同様に生成
            std::string gName = Form("GCheck_%s_Strip%02d", detNames[type], strip);
            if (hGainCheck_Global[gIdx]) delete hGainCheck_Global[gIdx];
            hGainCheck_Global[gIdx] = new TH1F(gName.c_str(), title.c_str(), 5000, 0, 100000); // こちらは粗くても良い
        }
    }
    
    // グラフ類の初期化
    hDeltaT_Neighbor = new TH1F("DeltaT_Nearest", "Delta T (Nearest);ns;Counts", 100000, 0, 100000);
    hDeltaT_n8 = new TH1F("DeltaT_n8", "Delta T (N to N+7);ns;Counts", 100000, 0, 100000);
    hGlobalPIndexMap = new TH1F("GlobalMap", "Hit Map;Global Index;Counts", 128, 0, 128);

    printLog("setupTree: Initialized 128 Histograms (1ns/bin resolution) in physical order.");
}

// ----------------------------------------------------------------------------
// 2. loadConversionFactors: マッピングのリンク (生成はしない)
// ----------------------------------------------------------------------------
void DetectorAnalyzer::loadConversionFactorsFromCSV(const std::string& csvFileName, const std::string& detName, int moduleID) {
    std::ifstream file(csvFileName);
    if (!file.is_open()) return;

    std::string line; 
    std::getline(file, line); // ヘッダースキップ
    
    activeModuleIDs_.insert(moduleID); 

    // このCSVファイルが担当する検出器タイプIDを特定
    int typeID = -1;
    if (detName == "X1") typeID = ID_X1;      
    else if (detName == "Y1") typeID = ID_Y1; 
    else if (detName == "X2") typeID = ID_X2; 
    else if (detName == "Y2") typeID = ID_Y2;
    
    if (typeID == -1) {
        printLog("[WARN] Unknown detector name in CSV load: " + detName);
        return;
    }

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line); 
        std::string cell; 
        std::vector<std::string> row;
        while (std::getline(ss, cell, ',')) row.push_back(cell);
        
        if (row.size() < 2) continue;

        try {
            int sysCh = std::stoi(row[0]);
            // Y検出器のチャンネルオフセット補正
            if (detName.find('Y') != std::string::npos) sysCh += Y_CHANNEL_OFFSET; 
            if (sysCh >= MAX_SYS_CH) continue;

            int localStrip = std::stoi(row[1]); // 1-32
            if (localStrip < 1 || localStrip > 32) continue;

            // LUT登録
            ch_LUT[moduleID][sysCh].isValid = true;
            ch_LUT[moduleID][sysCh].detTypeID = typeID;
            ch_LUT[moduleID][sysCh].strip = localStrip;
            
            // ★リンク処理: 先行生成されたGlobal配列から、該当するポインタを探してLUTに繋ぐ
            int gIdx = (typeID * 32) + (localStrip - 1);
            
            if (gIdx >= 0 && gIdx < 128) {
                // LUTにポインタをコピー（これでprocessChunkから高速アクセス可能）
                hRawToT_LUT[moduleID][sysCh] = hRawToT_Global[gIdx];
                hGainCheck_LUT[moduleID][sysCh] = hGainCheck_Global[gIdx];
                
                // Mapにも登録（保存やPDF生成で使う場合用）
                hRawToTMap[{moduleID, sysCh}] = hRawToT_Global[gIdx];
                gainCheckHists_[{moduleID, sysCh}] = hGainCheck_Global[gIdx];
                
                // コンフィグマップにも登録
                detConfigMap_[{moduleID, sysCh}] = {typeID, localStrip, detName, -1};
            }

        } catch (...) { continue; }
    }
    printLog("Mapped " + detName + " to Pre-allocated Histograms for Module " + std::to_string(moduleID));
}


// --- あかりさんのルールに基づくラン識別子抽出 ---
std::string DetectorAnalyzer::parseRunPrefix(const std::string& fullPath) {
    std::string fileName = std::filesystem::path(fullPath).filename().string();
    size_t firstUnderscore = fileName.find('_');
    if (firstUnderscore != std::string::npos) {
        return fileName.substr(0, firstUnderscore); // "PSD000001"
    }
    return fileName;
}

// --- フォルダ + PSD でランを特定する署名 ---
std::string DetectorAnalyzer::getRunSignature(const std::string& f) {
    std::filesystem::path p(f);
    return p.parent_path().string() + "###" + parseRunPrefix(f);
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

short DetectorAnalyzer::getRunID(const std::string& prefix) {
    if (runPrefixToId_.count(prefix)) return runPrefixToId_[prefix];
    short newId = (short)runIdToPrefix_.size();
    runIdToPrefix_.push_back(prefix); 
    runPrefixToId_[prefix] = newId;
    return newId;
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