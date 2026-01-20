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
#include "DetectorAnalyzer.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>

namespace fs = std::filesystem;

// ----------------------------------------------------------------------------
// グローバル変数 & ログ関数
// ----------------------------------------------------------------------------
static DetectorAnalyzer* g_currentAnalyzer = nullptr;


// ----------------------------------------------------------------------------
// グローバルログ関数（スパム防止版）
// ----------------------------------------------------------------------------
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
      tree_(nullptr), hDeltaT_Nearest(nullptr), hDeltaT_n8(nullptr), hGlobalPIndexMap(nullptr),
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
      currentRunDeathTime_(0),
      currentRunStartTime_(0),
      totalEffectiveTimeNs_(0),
      totalEffectiveTimeSec_(0.0),
      currentProcessingRunID_(-1),
      suppressTimeoutWarnings_(false), 
      prevEventTimeForEff_(0)
{
    g_currentAnalyzer = this; // グローバルポインタに登録
    gROOT->SetBatch(kTRUE); 

    // ★ 修正: Clock表示の起点をここで初期化
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


// ----------------------------------------------------------------------------
// デストラクタ
// ----------------------------------------------------------------------------
DetectorAnalyzer::~DetectorAnalyzer() {
    if (currentEventTime_ns_ > lastGlobalAnalysisTime_ns_) {
        analyzeGainShift(currentEventTime_ns_);
    }

    calculateEffectiveTime();

    if (outputFile_ && outputFile_->IsOpen()) {
        outputFile_->cd();
        if (tree_) tree_->Write();
        
        std::vector<std::pair<int, int>> sortedKeys;
        for (auto const& [key, info] : detConfigMap_) sortedKeys.push_back(key);
        std::sort(sortedKeys.begin(), sortedKeys.end(), [&](std::pair<int,int> a, std::pair<int,int> b){
            const auto& iA = detConfigMap_[a];
            const auto& iB = detConfigMap_[b];
            if(iA.detTypeID != iB.detTypeID) return iA.detTypeID < iB.detTypeID;
            return iA.strip < iB.strip; 
        });

        if (!outputFile_->GetDirectory("LogScalePlots")) outputFile_->mkdir("LogScalePlots");
        if (!outputFile_->GetDirectory("LinearScalePlots")) outputFile_->mkdir("LinearScalePlots");

        for (const auto& key : sortedKeys) {
            TH1F* hist = hRawToTMap[key];
            if (!hist || hist->GetEntries() < 10) continue;
            
            // 積分範囲の塗りつぶし用
            TH1F* hHL = nullptr;
            if (rangeMinMap_.count(key) && rangeMaxMap_.count(key)) {
                double minNS = (double)rangeMinMap_[key] * BIN_WIDTH_NS;
                double maxNS = (double)(rangeMaxMap_[key] + 1) * BIN_WIDTH_NS;
                
                int b1 = hist->FindBin(minNS);
                int b2 = hist->FindBin(maxNS);
                
                hHL = (TH1F*)hist->Clone(Form("%s_hl", hist->GetName()));
                hHL->Reset();
                for (int b = b1; b <= b2; ++b) hHL->SetBinContent(b, hist->GetBinContent(b));
                hHL->SetFillColor(kYellow); hHL->SetLineColor(kBlack);
            }

            outputFile_->cd("LinearScalePlots");
            TCanvas cLin(Form("c_lin_%s", hist->GetName()), hist->GetTitle());
            hist->Draw(); if (hHL) hHL->Draw("SAME HIST"); hist->Draw("SAME AXIS"); cLin.Write();

            outputFile_->cd("LogScalePlots");
            TCanvas cLog(Form("c_log_%s", hist->GetName()), hist->GetTitle());
            cLog.SetLogy(); hist->SetMinimum(0.5);
            hist->Draw(); if (hHL) hHL->Draw("SAME HIST"); hist->Draw("SAME AXIS"); cLog.Write();

            if (hHL) delete hHL;
        }

        outputFile_->cd();
        if(gT0Check) gT0Check->Write("T0_Check");

        if (!outputFile_->GetDirectory("GainHistory")) outputFile_->mkdir("GainHistory");
        outputFile_->cd("GainHistory");
        for (auto const& [key, g] : gainEvolutionGraphs_) if(g && g->GetN() > 0) g->Write();

        if (!outputFile_->GetDirectory("RateHistory")) outputFile_->mkdir("RateHistory");
        outputFile_->cd("RateHistory");
        for (auto const& [key, g] : rateEvolutionGraphs_) if(g && g->GetN() > 0) g->Write();

        outputFile_->Write();
    }
    
    if (outputFile_) { outputFile_->Close(); delete outputFile_; outputFile_ = nullptr; }
    if (gainLogCsv_.is_open()) gainLogCsv_.close();
    if (rateLogCsv_.is_open()) rateLogCsv_.close();
    printLog("Analysis cleanup completed.");
}

// ----------------------------------------------------------------------------
// データストリーム同期
// ----------------------------------------------------------------------------
bool DetectorAnalyzer::syncDataStream(std::ifstream& ifs, long long& skippedBytes) {
    skippedBytes = 0;
    
    // バッファサイズ (1MB)
    const size_t SCAN_BUFFER_SIZE = 1024 * 1024; 
    std::vector<char> buffer(SCAN_BUFFER_SIZE);
    
    long long globalOffset = ifs.tellg(); // 現在のファイル位置
    long long initialPos = globalOffset;
    
    // 最初のログだけ出す
    // printLog("[SYNC] Starting header search with buffered scan (1MB chunks)...");

    while (ifs) {
        ifs.read(buffer.data(), SCAN_BUFFER_SIZE);
        size_t bytesRead = ifs.gcount();
        if (bytesRead == 0) break;

        // バッファ内をスキャン
        for (size_t i = 0; i < bytesRead; ++i) {
            unsigned char h = static_cast<unsigned char>(buffer[i]);
            
            // ヘッダ候補が見つかった場合
            if (h == 0x69 || h == 0x6a) {
                bool syncOK = true;
                
                // バッファ内で確認できる場合
                if (i + 80 <= bytesRead) {
                    for (int k = 1; k < 10; ++k) {
                        unsigned char nextH = static_cast<unsigned char>(buffer[i + (k * 8)]);
                        if (nextH != 0x69 && nextH != 0x6a) {
                            syncOK = false;
                            break;
                        }
                    }
                } else {
                    // バッファ境界を跨ぐ場合: ファイルシークして慎重に確認
                    long long candidatePos = globalOffset + i;
                    long long nextChunkPos = globalOffset + bytesRead; // 読み込み後の位置
                    
                    for (int k = 1; k < 10; ++k) {
                         ifs.seekg(candidatePos + (k * 8), std::ios::beg);
                         char nextH;
                         if (!ifs.read(&nextH, 1)) { syncOK = false; break; }
                         if ((unsigned char)nextH != 0x69 && (unsigned char)nextH != 0x6a) {
                             syncOK = false; break;
                         }
                    }
                    
                    // ストリーム状態を復元（失敗時は次のチャンクへ、成功時はシークして終了）
                    if (!syncOK) {
                        ifs.clear();
                        ifs.seekg(nextChunkPos, std::ios::beg);
                    }
                }

                if (syncOK) {
                    // 同期成功！
                    long long finalPos = globalOffset + i;
                    ifs.seekg(finalPos, std::ios::beg);
                    skippedBytes = finalPos - initialPos;
                    // printLog("[SYNC] Sync found! Skipped bytes: " + std::to_string(skippedBytes));
                    return true;
                }
            }
        }

        // このチャンクで見つからなかった場合
        globalOffset += bytesRead;
        
        // 進捗ログ (10MBごとに表示)
        long long currentSkipped = globalOffset - initialPos;
        if (currentSkipped > 0 && currentSkipped % (10 * 1024 * 1024) == 0) {
             printLog("[SYNC] Scanning... Skipped: " + std::to_string(currentSkipped / 1024 / 1024) + " MB");
        }
    }
    
    // 全データ走査しても見つからなかった
    skippedBytes = globalOffset - initialPos;
    return false;
}


std::pair<unsigned long long, bool> DetectorAnalyzer::readEventsFromFile(const std::string& fileName, int modID, std::vector<Event>& rawEvents, long long& offset) {
    if (offset == 0) {
        printLog("[OPENING] Mod " + std::to_string(modID) + ": " + fileName);
        // currentRunPrefix_ は processBinaryFiles で管理するためここでは更新しない
    }

    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs.is_open()) {
        printLog("[ERROR] Could not open file: " + fileName);
        return {baseTimeMap_[modID], true};
    }

    ifs.seekg(0, std::ios::end);
    long long fileSize = ifs.tellg();
    ifs.seekg(offset, std::ios::beg);

    unsigned long long currentBaseTime = baseTimeMap_[modID]; 
    bool hasSeenHeader = hasSeenTimeHeaderMap_[modID]; 
    if (offset == 0) hasSeenHeader = false;

    const size_t BUFFER_SIZE = 64 * 1024 * 1024; 
    std::vector<char> buffer(BUFFER_SIZE);
    size_t leftover = 0;
    unsigned long long lastT = currentBaseTime; 

    ifs.read(buffer.data() + leftover, BUFFER_SIZE - leftover);
    size_t readCount = ifs.gcount();
    
    if (readCount > 0) {
        size_t totalBytesInBuffer = leftover + readCount;
        processedDataSize_ += readCount; 
        size_t i = 0;
        unsigned char* buf = reinterpret_cast<unsigned char*>(buffer.data());

        while (i + 8 <= totalBytesInBuffer) {
            unsigned char h = buf[i];
            bool syncOK = false;

            // 10個先までの同期チェック (既存ロジック維持)
            if (i + 80 > totalBytesInBuffer) {
                if (h == 0x69 || h == 0x6a) syncOK = true;
            } else {
                syncOK = true;
                for (int k = 1; k < 10; ++k) {
                    unsigned char nh = buf[i + (k * 8)];
                    if (nh != 0x69 && nh != 0x6a) { syncOK = false; break; }
                }
            }

            if (!hasSeenHeader) {
                if (syncOK && h == 0x69) {
                    unsigned char* p = &buf[i+1];
                    unsigned long long s = (static_cast<unsigned long long>(p[0]) << 24) | (p[1] << 16) | (p[2] << 8) | p[3];
                    unsigned long long ss = (static_cast<unsigned long long>(p[4]) << 2) | ((p[5] & 0xC0) >> 6);
                    unsigned long long t = ((p[5] & 0x3F) << 2) | ((p[6] & 0xC0) >> 6);
                    unsigned long long foundTime = s * 1000000000ULL + ss * 1000000ULL + t;

                    if (foundTime > baseTimeMap_[modID]) {
                        currentBaseTime = foundTime;
                        hasSeenHeader = true; 
                        lastT = currentBaseTime;
                        moduleAliveTime_[modID] = currentBaseTime;
                        if (firstT0Time_ == 0) firstT0Time_ = s;
                    }
                }
                i += (syncOK ? 8 : 1);
                continue; 
            }

            if (syncOK) {
                if (h == 0x69) { // T0 Header
                    unsigned char* p = &buf[i+1];
                    unsigned long long s = (static_cast<unsigned long long>(p[0]) << 24) | (p[1] << 16) | (p[2] << 8) | p[3];
                    unsigned long long ss = (static_cast<unsigned long long>(p[4]) << 2) | ((p[5] & 0xC0) >> 6);
                    unsigned long long t = ((p[5] & 0x3F) << 2) | ((p[6] & 0xC0) >> 6);
                    unsigned long long newTime = s * 1000000000ULL + ss * 1000000ULL + t;

                    long long diff = (long long)(newTime - currentBaseTime);
                    if (std::abs(diff) > 1000000000LL) {
                        if (diff > 0) { // Gap detected
                            deadTimeRanges_.push_back({currentBaseTime, newTime});
                            currentBaseTime = newTime;
                        }
                        // Glitch (diff < 0) is ignored
                    } else {
                        currentBaseTime = newTime;
                    }
                    lastT = currentBaseTime;
                    moduleAliveTime_[modID] = currentBaseTime;
                    i += 8;

                } else if (h == 0x6a) { // Event Data
                    unsigned char* p = &buf[i+1];
                    unsigned int tof = (p[0] << 16) | (p[1] << 8) | p[2];
                    unsigned int pw  = (p[3] << 12) | (p[4] << 4) | ((p[5] & 0xF0) >> 4);
                    int det = ((p[5] & 0x0F) << 8) | p[6];
                    int mod = (det >> 8) & 0xF;
                    int sys = det & 0xFF;
                    
                    if (mod < MAX_MODULES && sys < MAX_SYS_CH && ch_LUT[mod][sys].isValid) {
                        // ★物理補正 & 負値チェック
                        if (tof >= pw) { 
                            unsigned long long eventT = currentBaseTime + (unsigned long long)(tof - pw);
                            lastT = eventT;
                            moduleAliveTime_[modID] = eventT;
                            // ★RunID引数削除 (警告回避)
                            rawEvents.push_back({ ch_LUT[mod][sys].detTypeID, ch_LUT[mod][sys].strip, mod, eventT, (int)pw, sys });
                        }
                        // else: tof < pw なので捨てる (iだけ進める)
                    }
                    i += 8;
                } else { i += 8; }
            } else { i++; }
        }
        
        // ★頂いたコード通り、読み込み数分だけオフセットを進める
        offset += readCount; 
        
    } else { offset = fileSize; }

    baseTimeMap_[modID] = currentBaseTime;
    hasSeenTimeHeaderMap_[modID] = hasSeenHeader;
    return {lastT, (offset >= fileSize)}; 
}


void DetectorAnalyzer::processBinaryFiles(std::map<int, std::deque<std::string>>& fileQueues) {
    const size_t PER_MOD_CAP = 1000000; 
    std::map<int, std::deque<Event>> moduleBuffers;
    std::map<int, unsigned long long> moduleLastTime;
    auto lastUIDraw = std::chrono::system_clock::now();

    // 初期化
    for (auto const& [mod, q] : fileQueues) {
        moduleLastTime[mod] = 0; currentFileOffsets_[mod] = 0;
        hasSeenTimeHeaderMap_[mod] = false; baseTimeMap_[mod] = 0;
    }
    deadTimeRanges_.clear();
    isAnalysisStarted_ = false; globalStartTimestamp_ = 0;
    
    std::string currentRunSignature = ""; 
    unsigned long long lastSafeTime = 0;

    printSearchStatus();

    while (true) {
        // A. ラン開始判定 (バリア)
        if (currentRunSignature.empty()) {
            for (auto const& [m, q] : fileQueues) {
                if (!q.empty()) {
                    currentRunSignature = getRunSignature(q.front());
                    printLog("=== STARTING RUN: " + currentRunSignature + " ===");
                    currentRunStartTime_ = 0; lastGlobalAnalysisTime_ns_ = 0;
                    isAnalysisStarted_ = false; deadTimeRanges_.clear(); 
                    for(auto& [mod, t] : moduleLastTime) t = 0;
                    break; 
                }
            }
        }
        if (currentRunSignature.empty()) break;

        // B. 読み込み (ラン不一致なら待機)
        int laggingMod = -1;
        unsigned long long minLastTime = std::numeric_limits<unsigned long long>::max();
        for (auto const& [m, q] : fileQueues) {
            if (q.empty() || getRunSignature(q.front()) != currentRunSignature) continue;
            if (moduleLastTime[m] < minLastTime) { minLastTime = moduleLastTime[m]; laggingMod = m; }
        }

        if (laggingMod != -1 && moduleBuffers[laggingMod].size() < PER_MOD_CAP) {
            std::vector<Event> temp;
            auto res = readEventsFromFile(fileQueues[laggingMod].front(), laggingMod, temp, currentFileOffsets_[laggingMod]);
            for (auto& ev : temp) moduleBuffers[laggingMod].push_back(std::move(ev));
            moduleLastTime[laggingMod] = res.first;
            if (res.second) {
                printLog("[FINISHED] Mod " + std::to_string(laggingMod) + ": " + fileQueues[laggingMod].front());
                fileQueues[laggingMod].pop_front();
                currentFileOffsets_[laggingMod] = 0;
            }
        }

        // C. 同期・解析
        if (!isAnalysisStarted_) {
            bool ready = true; unsigned long long maxBase = 0;
            for (auto const& [m, q] : fileQueues) {
                if (!moduleBuffers[m].empty() || (!q.empty() && getRunSignature(q.front()) == currentRunSignature)) {
                    if (!hasSeenTimeHeaderMap_[m]) { ready = false; break; }
                    maxBase = std::max(maxBase, baseTimeMap_[m]);
                }
            }
            if (ready && maxBase > 0) {
                globalStartTimestamp_ = maxBase; lastSafeTime = maxBase;
                lastGlobalAnalysisTime_ns_ = maxBase; isAnalysisStarted_ = true;
                printLog("Sync Established.");
            }
        }

        if (isAnalysisStarted_) {
            unsigned long long safeTime = getSafeTime(moduleLastTime);
            if (safeTime > lastSafeTime) {
                if (safeTime - lastSafeTime < 10000000000ULL) totalEffectiveTimeNs_ += (safeTime - lastSafeTime);
                lastSafeTime = safeTime; currentEventTime_ns_ = safeTime;
                
                // ★デッドタイム情報の整理: 現在時刻より古い情報は不要なので捨てる (メモリ節約 & 高速化)
                while(!deadTimeRanges_.empty() && deadTimeRanges_.front().second < safeTime) {
                    deadTimeRanges_.erase(deadTimeRanges_.begin());
                }
            }
            std::vector<Event> mergeChunk;
            for (auto& [m, buf] : moduleBuffers) {
                while (!buf.empty() && buf.front().eventTime_ns <= safeTime) {
                    Event ev = std::move(buf.front()); buf.pop_front();
                    
                    // Veto処理
                    bool isDead = false;
                    for(const auto& r : deadTimeRanges_) {
                        if(ev.eventTime_ns >= r.first && ev.eventTime_ns <= r.second){ isDead=true; break;}
                    }
                    if (!isDead && ev.eventTime_ns >= globalStartTimestamp_) mergeChunk.push_back(std::move(ev));
                }
            }
            if (!mergeChunk.empty()) {
                std::sort(mergeChunk.begin(), mergeChunk.end(), [](const Event& a, const Event& b){ return a.eventTime_ns < b.eventTime_ns; });
                processChunk(mergeChunk);
            }
        }

        // D. バッファモニタ
        if (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - lastUIDraw).count() > 500) {
            printSearchStatus();
            std::stringstream ss; ss << "\n Run: " << currentRunSignature << "\n";
            for(auto const& [m, b] : moduleBuffers) {
                double usage = (double)b.size() / PER_MOD_CAP; if(usage>1.0) usage=1.0;
                ss << " Mod " << m << " [";
                for(int j=0; j<40; j++) ss << (j < usage*40 ? "█" : "░");
                ss << "] " << std::fixed << std::setprecision(1) << usage*100.0 << "%\n";
            }
            printLog(ss.str()); lastUIDraw = std::chrono::system_clock::now();
        }

        // E. ラン終了判定
        bool runOver = true;
        for (auto const& [m, q] : fileQueues) if (!q.empty() && getRunSignature(q.front()) == currentRunSignature) runOver = false;
        for (auto const& [m, b] : moduleBuffers) if (!b.empty()) runOver = false;
        if (runOver) {
            if(currentEventTime_ns_ > lastGlobalAnalysisTime_ns_) analyzeGainShift(currentEventTime_ns_);
            currentRunSignature = "";
        }
    }
}

bool DetectorAnalyzer::processChunk(const std::vector<Event>& sortedEvents) {
    if (sortedEvents.empty()) return true;

    for (size_t i = 0; i < sortedEvents.size(); ++i) {
        const auto& e = sortedEvents[i];
        
        // --- 削除するブロック (ここが runID を使っていた) ---
        // if (currentProcessingRunID_ != -1 && e.runID != currentProcessingRunID_) { ... }
        // currentProcessingRunID_ = e.runID;
        // ------------------------------------------------

        // ラン開始時の初期化 (ここは runID に依存せず tot で判定しているので維持)
        if (currentRunStartTime_ == 0) {
            if (e.tot >= 1000) {
                currentRunStartTime_ = e.eventTime_ns;
                lastTimeoutCheckTime_ = e.eventTime_ns; 
                prevEventTimeForEff_ = e.eventTime_ns;
                for(int mid : activeModuleIDs_) moduleAliveTime_[mid] = e.eventTime_ns;
            } else continue;
        }

        currentEventTime_ns_ = e.eventTime_ns;
        moduleAliveTime_[e.module] = e.eventTime_ns;

        // システム健全性チェックと実効時間加算
        bool isSystemHealthy = true;
        for (int mid : activeModuleIDs_) {
            if (moduleAliveTime_.find(mid) == moduleAliveTime_.end()) {
                isSystemHealthy = false; break;
            }
        }
        if (isSystemHealthy && prevEventTimeForEff_ > 0) {
            unsigned long long dt_ns = e.eventTime_ns - prevEventTimeForEff_;
            if (dt_ns < 5000000000ULL) totalEffectiveTimeNs_ += dt_ns;
        }
        prevEventTimeForEff_ = e.eventTime_ns;

        // ゲイン解析トリガー
        if (e.eventTime_ns >= lastGlobalAnalysisTime_ns_ + GAIN_ANALYSIS_WINDOW_NS) {
            analyzeGainShift(e.eventTime_ns);
            lastGlobalAnalysisTime_ns_ = e.eventTime_ns;
        }

        // データ充填
        foundChannels_[{e.module, e.sysCh}].insert(e.strip);
        if (e.module < MAX_MODULES && e.sysCh < MAX_SYS_CH) {
            if (hGainCheck_LUT[e.module][e.sysCh]) hGainCheck_LUT[e.module][e.sysCh]->Fill((double)e.tot);
            if (hRawToT_LUT[e.module][e.sysCh]) hRawToT_LUT[e.module][e.sysCh]->Fill((double)e.tot);
        }

        b_type = e.type; b_strip = e.strip; b_time = (long long)e.eventTime_ns; b_tot = e.tot;
        if (tree_) tree_->Fill();

        // N+7 相関解析
        if (i > 0) {
            long long dt1 = (long long)e.eventTime_ns - (long long)sortedEvents[i-1].eventTime_ns;
            if (dt1 >= 0 && dt1 < 1000000) hDeltaT_Nearest->Fill((double)dt1);
        }
        if (i + 7 < sortedEvents.size()) {
            long long dt8 = (long long)sortedEvents[i+7].eventTime_ns - (long long)e.eventTime_ns;
            if (dt8 >= 0 && dt8 < 1000000) hDeltaT_n8->Fill((double)dt8);
        }
    }
    return true;
}


void DetectorAnalyzer::printSearchStatus() {
    static size_t lastReportedCount = 0;

    if (!isAnalysisStarted_) {
        std::cout << "\r\033[K" << "[STATUS] Syncing..." << std::flush;
        return;
    }

    auto now = std::chrono::system_clock::now();
    double elapsedClock = std::chrono::duration_cast<std::chrono::seconds>(now - analysisStartTime_).count();
    
    // LiveTime と DataTime は同じ変数(totalEffectiveTimeNs_)を使う
    double liveTimeMin = (double)totalEffectiveTimeNs_ / 1.0e9 / 60.0;
    double progress = (analysisDuration_ns_ > 0) ? (double)totalEffectiveTimeNs_ / (double)analysisDuration_ns_ * 100.0 : 0.0;
    
    // ★重要修正: 100%キャップを撤廃 (事実をそのまま表示)

    std::cout << "\r\033[K" 
              << "Progress: [" << std::fixed << std::setprecision(1) << progress << "%] "
              << "Data: " << liveTimeMin << " min | " 
              << "Real: " << (int)elapsedClock / 60 << "m" << (int)elapsedClock % 60 << "s | "
              << "Live: " << liveTimeMin << " min " // 単位を min に統一
              << std::flush;

    size_t currentCount = 0;
    for (auto const& [key, strips] : foundChannels_) currentCount += strips.size();

    if (currentCount > lastReportedCount) {
        std::cout << "\n\n[NEW CHANNELS DETECTED] Total: " << currentCount << " (+ " << (currentCount - lastReportedCount) << ")" << std::endl;
        
        std::map<int, std::set<int>> activeStripsByDet;
        for (auto const& [key, strips] : foundChannels_) {
            int mod = key.first; int sys = key.second;
            if (ch_LUT[mod][sys].isValid) {
                for(int s : strips) activeStripsByDet[ch_LUT[mod][sys].detTypeID].insert(s);
            }
        }
        std::string detNames[] = {"X1", "Y1", "X2", "Y2"};
        for (int d = 0; d < 4; ++d) {
            std::cout << "  " << detNames[d] << ": ";
            const auto& activeSet = activeStripsByDet[d];
            for (int s = 1; s <= 32; ++s) {
                if (activeSet.count(s)) std::cout << "\033[32m" << std::setw(2) << s << "\033[0m ";
                else std::cout << "\033[90m..\033[0m ";
            }
            std::cout << std::endl;
        }
        std::cout << "--------------------------------------------------------------------------------" << std::endl;
        lastReportedCount = currentCount;
    }
}


// 解析停止直前のバイナリダンプ
// ----------------------------------------------------------------------------
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

unsigned long long DetectorAnalyzer::getSafeTime(const std::map<int, unsigned long long>& lastTimes) {
    unsigned long long minT = std::numeric_limits<unsigned long long>::max();
    // lastTimesが空の場合を考慮
    if (lastTimes.empty()) return 0;
    for (auto const& [id, t] : lastTimes) {
        if (t < minT) minT = t;
    }
    return minT;
}



// ----------------------------------------------------------------------------
// ゲイン解析ロジック
// ----------------------------------------------------------------------------
int DetectorAnalyzer::findRightMostPeak(TH1F* hist) {
    if (!hist) return 0;
    hist->GetXaxis()->SetRangeUser(PEAK_SEARCH_MIN_TOT, PEAK_SEARCH_MAX_TOT);
    int maxBin = hist->GetMaximumBin();
    hist->GetXaxis()->SetRange(0, 0); 
    return maxBin;
}


double DetectorAnalyzer::calculateResidual(TH1F* hTarget, TH1F* hTemplate, int shiftBins, int binMin, int binMax) {
    int tBinMin = binMin + shiftBins;
    int tBinMax = binMax + shiftBins;
    if (tBinMin < 1 || tBinMax > hTarget->GetNbinsX()) return -1.0;
    
    double intTarget = hTarget->Integral(tBinMin, tBinMax);
    double intTemplate = hTemplate->Integral(binMin, binMax);
    if (intTarget <= 0 || intTemplate <= 0) return -1.0; 
    
    double scale = intTemplate / intTarget;
    double residualSum = 0.0;
    for (int i = binMin; i <= binMax; ++i) {
        int targetBin = i + shiftBins; 
        double diff = (hTarget->GetBinContent(targetBin) * scale) - hTemplate->GetBinContent(i);
        residualSum += (diff * diff);
    }
    return residualSum;
}


// 5us (250bin) 空白探索
void DetectorAnalyzer::determineIntegrationRange(TH1F* hist, int peakBin, int& outMinBin, int& outMaxBin) {
    int nBins = hist->GetNbinsX();
    const int SILENCE_BINS = 250; 
    
    outMinBin = 1; int zeroCount = 0;
    for (int i = peakBin; i > 1; --i) {
        if (hist->GetBinContent(i) <= 0) zeroCount++; else zeroCount = 0;
        if (zeroCount >= SILENCE_BINS) { outMinBin = i + zeroCount; if(outMinBin<1)outMinBin=1; break; }
    }
    outMaxBin = nBins; zeroCount = 0;
    for (int i = peakBin; i < nBins; ++i) {
        if (hist->GetBinContent(i) <= 0) zeroCount++; else zeroCount = 0;
        if (zeroCount >= SILENCE_BINS) { outMaxBin = i - zeroCount; if(outMaxBin>nBins)outMaxBin=nBins; break; }
    }
}

// マッチング (端に当たったらエラー)
int DetectorAnalyzer::findBestShift(TH1F* hTarget, TH1F* hTemplate, int binMin, int binMax, const std::string&) {
    double minRes = 1e20; int bestS = SHIFT_CALC_ERROR; bool found = false;
    const int SEARCH_RANGE = 25; // ±500ns

    for (int s = -SEARCH_RANGE; s <= SEARCH_RANGE; ++s) { 
        double res = calculateResidual(hTarget, hTemplate, s, binMin, binMax);
        if (res >= 0 && res < minRes) { minRes = res; bestS = s; found = true; }
    }
    if (!found || std::abs(bestS) == SEARCH_RANGE) return SHIFT_CALC_ERROR;
    return bestS;
}

bool DetectorAnalyzer::analyzeGainShift(unsigned long long currentTime_ns) { 
    if (!gainLogCsv_.is_open()) return true;
    if (!isCsvHeaderWritten_) {
        // (ヘッダ書き込み: 既存ママ)
        struct SortInfo { std::pair<int,int> k; int t; int s; std::string n; };
        std::vector<SortInfo> sorted;
        for (auto const& [keyPair, info] : detConfigMap_) sorted.push_back({keyPair, info.detTypeID, info.strip, info.cachedName});
        std::sort(sorted.begin(), sorted.end(), [](const SortInfo& a, const SortInfo& b){ if(a.t != b.t) return a.t < b.t; return a.s < b.s; });
        gainLogCsv_ << "Time_ns"; rateLogCsv_ << "Time_ns";
        for(const auto& s : sorted) { std::string label = s.n + "-" + std::to_string(s.s); gainLogCsv_ << "," << label; rateLogCsv_ << "," << label; }
        gainLogCsv_ << "\n"; rateLogCsv_ << "\n";
        isCsvHeaderWritten_ = true;
    }

    gainLogCsv_ << currentTime_ns; rateLogCsv_ << currentTime_ns;
    double durMin = (double)GAIN_ANALYSIS_WINDOW_NS / 1.0e9 / 60.0; // 10分

    std::vector<std::pair<int, int>> sortedKeys;
    for (auto const& [key, info] : detConfigMap_) sortedKeys.push_back(key);
    std::sort(sortedKeys.begin(), sortedKeys.end(), [&](std::pair<int,int> a, std::pair<int,int> b){
        const auto& iA = detConfigMap_[a]; const auto& iB = detConfigMap_[b];
        if(iA.detTypeID != iB.detTypeID) return iA.detTypeID < iB.detTypeID; return iA.strip < iB.strip; 
    });

    for (const auto& key : sortedKeys) {
        TH1F* hist = hGainCheck_LUT[key.first][key.second];
        std::string sVal = "ND"; double rVal = 0.0;

        if (hist) {
            // カウントによるND判定は行わない (ご指示通り)
            
            if (!isTemplateCaptured_) {
                // 初回
                int pBin = findRightMostPeak(hist); int rMin, rMax;
                determineIntegrationRange(hist, pBin, rMin, rMax);
                rangeMinMap_[key] = rMin; rangeMaxMap_[key] = rMax;
                
                templateHists_[key] = (TH1F*)hist->Clone(Form("%s_tmpl", hist->GetName()));
                templateHists_[key]->SetDirectory(0);
                
                cumulativeShiftMap_[key] = 0.0; sVal = "0";
                if(rMin>0 && rMax>rMin) rVal = hist->Integral(rMin, rMax) / durMin;
            } else {
                // 2回目以降
                int dBin = findBestShift(hist, templateHists_[key], rangeMinMap_[key], rangeMaxMap_[key], "");
                
                if (dBin != SHIFT_CALC_ERROR) {
                    if (dBin != 0) {
                        rangeMinMap_[key] += dBin; rangeMaxMap_[key] += dBin;
                        cumulativeShiftMap_[key] += (double)dBin * 20.0;
                        delete templateHists_[key];
                        templateHists_[key] = (TH1F*)hist->Clone(Form("%s_tmpl", hist->GetName()));
                        templateHists_[key]->SetDirectory(0);
                    }
                    sVal = std::to_string(cumulativeShiftMap_[key]);
                    rVal = hist->Integral(rangeMinMap_[key], rangeMaxMap_[key]) / durMin;
                } else {
                    // マッチング失敗時のみ ND
                    sVal = "ND"; 
                    rVal = 0.0; // あるいは前回値を維持するか、NDなら0にするか。ここでは0にします。
                }
            }

            if(sVal != "ND") {
                if(gainEvolutionGraphs_.count(key)) gainEvolutionGraphs_[key]->SetPoint(gainEvolutionGraphs_[key]->GetN(), (double)currentTime_ns, cumulativeShiftMap_[key]);
                if(rateEvolutionGraphs_.count(key)) rateEvolutionGraphs_[key]->SetPoint(rateEvolutionGraphs_[key]->GetN(), (double)currentTime_ns, rVal);
            }
            hist->Reset();
        }
        gainLogCsv_ << "," << sVal; rateLogCsv_ << "," << rVal;
    }
    gainLogCsv_ << std::endl; rateLogCsv_ << std::endl;
    if(!isTemplateCaptured_) isTemplateCaptured_ = true;
    return true; 
}


// ----------------------------------------------------------------------------
// ユーティリティ
// ----------------------------------------------------------------------------

std::string DetectorAnalyzer::getRunSignature(const std::string& fileName) {
    namespace fs = std::filesystem;
    fs::path p(fileName);
    std::string parentDir = p.parent_path().string();
    std::string filenameStr = p.filename().string();
    
    static std::regex run_re(R"((PSD\d+_\d+))");
    std::smatch match;
    std::string fileRunID = "UnknownID";
    if (std::regex_search(filenameStr, match, run_re)) {
        fileRunID = match[1].str();
    }
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
    struct FileInfo {
        long long date; 
        int runID; 
        int modID;
        int fileSeq; 
        std::string path;
        
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

        long long fDate = 0; 
        int fRun = 0;
        
        size_t pSep = parentDir.find('_');
        if (pSep != std::string::npos) {
            try {
                if (parentDir.substr(0, 3) == "PSD") {
                    fRun = std::stoi(parentDir.substr(3, pSep - 3));
                }
                fDate = std::stoll(parentDir.substr(pSep + 1));
            } catch (...) {}
        }

        std::vector<std::string> tokens;
        std::stringstream ss(fileName); 
        std::string item;
        while (std::getline(ss, item, '_')) { tokens.push_back(item); }

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

    // キューへの詰め込み
    for (const auto& info : allFiles) {
        fileQueues[info.modID].push_back(info.path);
    }

    // --- 復活させた詳細ログ出力 ---
    printLog("================================================================");
    printLog("[INFO] Final Sorted File Queue Listing:");
    for (auto const& [mod, q] : fileQueues) {
        printLog("Module " + std::to_string(mod) + ": " + std::to_string(q.size()) + " files found.");
        for (const auto& path : q) {
            printLog("  [Queue Mod " + std::to_string(mod) + "] " + path);
        }
    }
    printLog("================================================================");
}


void DetectorAnalyzer::loadConversionFactorsFromCSV(const std::string& csvFileName, const std::string& detName, int moduleID) {
    std::ifstream file(csvFileName);
    if (!file.is_open()) return;

    int typeID = -1;
    if (detName == "X1") typeID = ID_X1; else if (detName == "Y1") typeID = ID_Y1;
    else if (detName == "X2") typeID = ID_X2; else if (detName == "Y2") typeID = ID_Y2;

    std::string line;
    std::getline(file, line); 
    
    if (nameStripToMapKey_.find(detName) == nameStripToMapKey_.end()) {
        nameStripToMapKey_[detName]; 
    }
    
    activeModuleIDs_.insert(moduleID); 

    struct CSVItem { int localCh; std::pair<int, int> key; std::string label; };
    std::vector<CSVItem> sortedItems;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;
        while (std::getline(ss, cell, ',')) { row.push_back(cell); }
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

            if (offlineStrips_.count({typeID, localCh})) {
                ignoredMonitorChannels_.insert(keyPair);
            }

            if (moduleID < MAX_MODULES) {
                ch_LUT[moduleID][sysCh].isValid = true;
                ch_LUT[moduleID][sysCh].detTypeID = typeID;
                ch_LUT[moduleID][sysCh].strip = localCh;
                ch_LUT[moduleID][sysCh].pIndex = pIndex;
            }

        } catch (...) { }
    }

    std::sort(sortedItems.begin(), sortedItems.end(), [](const CSVItem& a, const CSVItem& b) {
        return a.localCh < b.localCh;
    });

    outputFile_->cd(); 
    for (const auto& item : sortedItems) {
        std::string hName = "RawToT_" + item.label;
        std::string hTitle = "Raw ToT: " + item.label;
        if (gDirectory->Get(hName.c_str())) delete gDirectory->Get(hName.c_str());
        
        if (hRawToTMap.find(item.key) == hRawToTMap.end()) {
            hRawToTMap[item.key] = new TH1F(hName.c_str(), hTitle.c_str(), 100000, 0, 100000);
        }
        
        if (gainCheckHists_.find(item.key) == gainCheckHists_.end()) {
            std::string gName = Form("GCheck_%s", item.label.c_str());
            if(gDirectory->Get(gName.c_str())) delete gDirectory->Get(gName.c_str());
            gainCheckHists_[item.key] = new TH1F(gName.c_str(), Form("Gain Check %s", item.label.c_str()), MONITOR_HIST_BINS, 0, MONITOR_HIST_MAX_TOT); 
            gainCheckHists_[item.key]->SetDirectory(0); 
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

    std::string line;
    std::getline(file, line); 
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;
        while (std::getline(ss, cell, ',')) { row.push_back(cell); }

        if (row.size() < 3) continue;
        try {
            int detType = std::stoi(row[0]);
            int stripNo = std::stoi(row[1]);
            int flag    = std::stoi(row[2]);

            if (flag == 1) {
                offlineStrips_.insert({detType, stripNo});
            }
        } catch (...) {}
    }
}

void DetectorAnalyzer::setTimeLimitMinutes(double minutes) {
    if (minutes <= 0) {
        analysisDuration_ns_ = std::numeric_limits<unsigned long long>::max();
        printLog("Time limit set to: UNLIMITED");
    } else {
        analysisDuration_ns_ = (unsigned long long)(minutes * 60.0 * 1000000000.0);
        printLog("Time limit set to: " + std::to_string(minutes) + " minutes");
    }
}

void DetectorAnalyzer::setupTree() {
    outputFile_->cd();
    tree_ = new TTree("Events", "Raw Detector Hits");
    tree_->Branch("type",  &b_type,  "type/I");   
    tree_->Branch("strip", &b_strip, "strip/I");  
    tree_->Branch("time",  &b_time,  "time/L");   
    tree_->Branch("tot",   &b_tot,   "tot/I");    

    hDeltaT_Nearest = new TH1F("DeltaT_Nearest", "Delta T (Nearest Neighbor);Delta T [ns];Counts / ns", 100000, 0, 100000);
    hDeltaT_n8 = new TH1F("DeltaT_n_plus_7", "Delta T (N to N+7);Delta T [ns];Counts / ns", 100000, 0, 100000);
    hGlobalPIndexMap = new TH1F("GlobalHitMap_PIndex", "Global Hit Count vs P-Index;P-index;Total Counts", 600, 0, 600);
    
    gT0Check = new TGraph();
    gT0Check->SetName("T0_Check");
    gT0Check->SetTitle("T0 Linearity Check;Count;Time Elapsed [s]");
    gT0Check->SetMarkerStyle(20);
    gT0Check->SetMarkerSize(0.5);
    gT0Check->SetMarkerColor(kRed);
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