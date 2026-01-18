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

namespace fs = std::filesystem;

// ----------------------------------------------------------------------------
// グローバル変数 & ログ関数
// ----------------------------------------------------------------------------
static DetectorAnalyzer* g_currentAnalyzer = nullptr;

void printLog(const std::string& msg) {
    // 画面クリア (スクロールさせるため、現在行をクリアして上書き風に見せる)
    std::cout << "\r\033[K" << std::flush;

    // ログ出力
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm* now_tm = std::localtime(&now_c);
    std::cout << "[" << std::put_time(now_tm, "%Y-%m-%d %H:%M:%S") << "] [INFO] " << msg << std::endl;

    // ダッシュボードを再描画して、進捗が消えないようにする
    if (g_currentAnalyzer) {
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

// ----------------------------------------------------------------------------
// ファイル読み込み (修正済み: EOF無限ループ防止 + ハイブリッドスキップ)
// ----------------------------------------------------------------------------
std::pair<unsigned long long, bool> DetectorAnalyzer::readEventsFromFile(
    const std::string& fileName, int modID, std::vector<Event>& rawEvents, long long& offset
) {
    if (offset == 0) {
        printLog("Opening NEW file [Mod " + std::to_string(modID) + "]: " + fileName);
    }
    
    std::string prefix = getRunSignature(fileName);
    short rid = getRunID(prefix); 
    currentRunPrefix_[modID] = prefix;

    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs.is_open()) return {0, true}; 

    unsigned long long currentBaseTime = baseTimeMap_[modID]; 
    bool hasSeenHeader = hasSeenTimeHeaderMap_[modID]; 

    // --- オフセット適用と初期Sync ---
    if (offset == 0) {
        long long skip = 0;
        if (!syncDataStream(ifs, skip)) {
            printLog("[WARNING] Mod " + std::to_string(modID) + ": No sync found in header area.");
            return {0, true};
        }
        offset = ifs.tellg(); 
        hasSeenHeader = false; 
    } else {
        ifs.seekg(offset, std::ios::beg);
    }

    ifs.seekg(0, std::ios::end);
    long long fileSize = ifs.tellg();
    ifs.seekg(offset, std::ios::beg);

    const size_t BUFFER_SIZE = 64 * 1024 * 1024; 
    std::vector<char> buffer(BUFFER_SIZE);
    const size_t MAX_EVENTS_IN_THIS_CALL = 2000000; 
    size_t eventsReadThisCall = 0;
    size_t leftover = 0;
    
    unsigned long long lastT = currentBaseTime; 
    bool timeLimitReached = false;

    while (ifs && eventsReadThisCall < MAX_EVENTS_IN_THIS_CALL) {
        ifs.read(buffer.data() + leftover, BUFFER_SIZE - leftover);
        size_t readCount = ifs.gcount();
        
        // ★★★ 修正の核心: 無限ループ防止 ★★★
        // 読み込みサイズが0で、残りが1パケット(8バイト)未満なら、これ以上読めないため強制終了
        if (readCount == 0 && leftover < 8) {
            offset = fileSize; // ファイルサイズまで進めて「完了」扱いにし、ループを脱出
            break; 
        }

        size_t totalBytesInBuffer = leftover + readCount;
        size_t i = 0;
        unsigned char* buf = reinterpret_cast<unsigned char*>(buffer.data());

        while (i + 8 <= totalBytesInBuffer) {
            unsigned char h = buf[i];

            // モードA: T0探索 (ハイブリッド・スキップ)
            if (!hasSeenHeader) {
                if (h == 0x69) { // T0候補
                    if (i + 80 > totalBytesInBuffer) break; 
                    bool ok = true;
                    for (int k = 1; k < 10; ++k) {
                        unsigned char nh = buf[i + (k * 8)];
                        if (nh != 0x69 && nh != 0x6a) { ok = false; break; }
                    }
                    if (ok) {
                        unsigned char* p = &buf[i+1];
                        unsigned long long s = (static_cast<unsigned long long>(p[0]) << 24) | (p[1] << 16) | (p[2] << 8) | p[3];
                        unsigned long long ss = (static_cast<unsigned long long>(p[4]) << 2) | ((p[5] & 0xC0) >> 6);
                        unsigned long long t = ((p[5] & 0x3F) << 2) | ((p[6] & 0xC0) >> 6);
                        
                        if (ss < 1000 && t < 50) {
                            currentBaseTime = s * 1000000000ULL + ss * 1000000ULL + t;
                            hasSeenHeader = true;
                            lastT = currentBaseTime;
                            moduleAliveTime_[modID] = currentBaseTime;
                            i += 8; continue;
                        }
                    }
                } 
                else if (h == 0x6a) { 
                    i += 8; continue; // データパケット高速スキップ
                }
                i++; continue; // ゴミなら1バイト
            }

            // モードB: 通常解析
            bool syncOK = true;
            if (i + 80 <= totalBytesInBuffer) {
                for (int k = 1; k < 10; ++k) {
                    unsigned char nh = buf[i + (k * 8)];
                    if (nh != 0x69 && nh != 0x6a) { syncOK = false; break; }
                }
            }

            if (syncOK) {
                if (h == 0x69) {
                    unsigned char* p = &buf[i+1];
                    unsigned long long s = (static_cast<unsigned long long>(p[0]) << 24) | (p[1] << 16) | (p[2] << 8) | p[3];
                    unsigned long long ss = (static_cast<unsigned long long>(p[4]) << 2) | ((p[5] & 0xC0) >> 6);
                    unsigned long long t = ((p[5] & 0x3F) << 2) | ((p[6] & 0xC0) >> 6);
                    if (ss < 1000 && t < 50) {
                        unsigned long long newTime = s * 1000000000ULL + ss * 1000000ULL + t;
                        if (newTime < currentBaseTime + 100000000000ULL) {
                            currentBaseTime = newTime; 
                            lastT = currentBaseTime;
                            moduleAliveTime_[modID] = currentBaseTime;
                        }
                    }
                    i += 8;
                } else if (h == 0x6a) {
                    unsigned char* p = &buf[i+1];
                    unsigned int tof = (p[0] << 16) | (p[1] << 8) | p[2];
                    unsigned int pw  = (p[3] << 12) | (p[4] << 4) | ((p[5] & 0xF0) >> 4);
                    int det = ((p[5] & 0x0F) << 8) | p[6];
                    int mod = (det >> 8) & 0xF;
                    int sys = det & 0xFF;

                    if (mod < MAX_MODULES && sys < MAX_SYS_CH && ch_LUT[mod][sys].isValid) {
                        long long diff = (long long)tof - (long long)pw;
                        if (!(diff < 0 && currentBaseTime < (unsigned long long)std::abs(diff))) {
                            unsigned long long eventT = currentBaseTime + diff;
                            lastT = eventT;
                            moduleAliveTime_[modID] = eventT;

                            if (globalStartTimestamp_ == 0 || eventT < globalStartTimestamp_) globalStartTimestamp_ = eventT;
                            if (analysisDuration_ns_ != std::numeric_limits<unsigned long long>::max() && eventT > globalStartTimestamp_ + analysisDuration_ns_) {
                                timeLimitReached = true;
                            } else {
                                rawEvents.push_back({ ch_LUT[mod][sys].detTypeID, ch_LUT[mod][sys].strip, mod, eventT, (int)pw, sys, rid });
                                foundChannels_[{mod, sys}].insert(ch_LUT[mod][sys].strip); 
                                eventsReadThisCall++;
                            }
                        }
                    }
                    i += 8;
                } else {
                    i += 8;
                }
            } else {
                i++; // 同期喪失時は1バイトシフト
            }
        }
        
        offset += i;
        processedDataSize_ += i; 

        if (timeLimitReached) break; 
        
        leftover = totalBytesInBuffer - i;
        if (leftover > 0) std::memmove(buffer.data(), buffer.data() + i, leftover);
        if (offset >= fileSize) break;
    }
    
    baseTimeMap_[modID] = currentBaseTime;
    hasSeenTimeHeaderMap_[modID] = hasSeenHeader;
    return {lastT, (offset >= fileSize || timeLimitReached)}; 
}

// ----------------------------------------------------------------------------
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
    for (auto const& [id, t] : lastTimes) { if (t < minT) minT = t; }
    return minT;
}

// ----------------------------------------------------------------------------
// メイン解析ループ
// ----------------------------------------------------------------------------

void DetectorAnalyzer::processBinaryFiles(std::map<int, std::deque<std::string>>& fileQueues) {
    const size_t BUFFER_CAP = 2000000;
    std::vector<Event> eventChunk;
    eventChunk.reserve(BUFFER_CAP * 2);
    
    std::map<int, unsigned long long> moduleLastTime;
    
    // 初期化
    for (auto const& [mod, q] : fileQueues) {
        moduleLastTime[mod] = 0; 
        currentFileOffsets_[mod] = 0;
        baseTimeMap_[mod] = 0; 
        hasSeenTimeHeaderMap_[mod] = false;
        moduleAliveTime_[mod] = 0;
        
        if (!q.empty()) {
            currentRunPrefix_[mod] = getRunSignature(q.front());
            printLog("[INFO] Mod " + std::to_string(mod) + ": Initial file -> " + q.front());
        }
    }

    printLog("Step 3: Parallel streaming analysis started.");
    analysisStartTime_ = std::chrono::system_clock::now();
    lastGlobalAnalysisTime_ns_ = 0; 

    // デバッグ: ストール（膠着）検知用カウンタ
    int stallCounter = 0;

    while (true) {
        bool allFilesFinished = true;
        bool dataReadInThisLoop = false;
        bool dataProcessedInThisLoop = false;

        bool anyQueueRemaining = false;
        for (auto const& [modID, queue] : fileQueues) {
            if (!queue.empty()) anyQueueRemaining = true;
        }
        if (!anyQueueRemaining) {
            printLog("[DEBUG] Stopping analysis: All queues are empty.");
            break; 
        }

        // SafeTime (全モジュールの最小時刻) の計算
        unsigned long long currentMinTime = getSafeTime(moduleLastTime);

        for (auto& [modID, queue] : fileQueues) {
            if (!queue.empty()) allFilesFinished = false;
            
            // -------------------------------------------------------------
            // スケジューラ判定デバッグ
            // -------------------------------------------------------------
            bool bufferOK = eventChunk.size() < BUFFER_CAP;
            bool timeOK   = moduleLastTime[modID] <= currentMinTime;
            bool shouldRead = !queue.empty() && (bufferOK || timeOK);

            /*
            // 膠着時のみ出す詳細デバッグ（普段はコメントアウト推奨だが、今回はONにするか検討）
            if (stallCounter > 5) {
                std::string reason = "";
                if (queue.empty()) reason += "QueueEmpty ";
                if (!bufferOK) reason += "BufferFull ";
                if (!timeOK) reason += "TimeAhead(" + std::to_string(moduleLastTime[modID]) + " > " + std::to_string(currentMinTime) + ") ";
                printLog("[DEBUG_STALL] Mod " + std::to_string(modID) + " Skip Reason: " + reason);
            }
            */

            if (shouldRead) {
                // 読み込み実行
                auto res = readEventsFromFile(queue.front(), modID, eventChunk, currentFileOffsets_[modID]);
                
                if (res.first > 0) moduleLastTime[modID] = res.first;
                dataReadInThisLoop = true;
                
                // ファイル完了処理
                if (res.second) { 
                    printLog("[INFO] Mod " + std::to_string(modID) + ": Finished reading " + queue.front());
                    
                    queue.pop_front(); 
                    currentFileOffsets_[modID] = 0; 
                    
                    if (!queue.empty()) {
                        printLog("[INFO] Mod " + std::to_string(modID) + ": Next file -> " + queue.front());
                    } else {
                        printLog("[INFO] Mod " + std::to_string(modID) + ": Queue empty. No more files.");
                    }

                    if (analysisDuration_ns_ != std::numeric_limits<unsigned long long>::max() && 
                        moduleLastTime[modID] > globalStartTimestamp_ + analysisDuration_ns_) {
                         printLog("[INFO] Global time limit reached via Module " + std::to_string(modID));
                         allFilesFinished = true; 
                         goto END_ANALYSIS; 
                    }
                }
            }
        }

        // イベントソートと処理
        if (!eventChunk.empty()) {
            std::sort(eventChunk.begin(), eventChunk.end(), [](const Event& a, const Event& b){
                return a.eventTime_ns < b.eventTime_ns;
            });

            unsigned long long safeTime = std::numeric_limits<unsigned long long>::max();
            bool hasActive = false;
            for (auto const& [m, q] : fileQueues) {
                if (!q.empty()) { safeTime = std::min(safeTime, moduleLastTime[m]); hasActive = true; }
            }
            if (!hasActive || allFilesFinished) safeTime = eventChunk.back().eventTime_ns + 1;

            auto it = std::upper_bound(eventChunk.begin(), eventChunk.end(), safeTime, 
                [](unsigned long long t, const Event& e){ return t < e.eventTime_ns; });
            size_t cut = std::distance(eventChunk.begin(), it);

            if (cut > 0) {
                std::vector<Event> safe(eventChunk.begin(), eventChunk.begin() + cut);
                processChunk(safe);
                eventChunk.erase(eventChunk.begin(), eventChunk.begin() + cut);
                dataProcessedInThisLoop = true;
            }
        }
        
        if (allFilesFinished && eventChunk.empty()) break;
        
        // -------------------------------------------------------------
        // デッドロック（空回り）検知 & 強制ダンプ
        // -------------------------------------------------------------
        if (!dataReadInThisLoop && !dataProcessedInThisLoop && !allFilesFinished) {
            stallCounter++;
            
            // 1万回ループが空回りしたら（=フリーズ）、状況をダンプする
            if (stallCounter % 10000 == 0) {
                printLog("!!! STALL DETECTED !!! Loop cycling without work.");
                printLog("Global SafeTime: " + std::to_string(currentMinTime));
                printLog("Buffer Size: " + std::to_string(eventChunk.size()) + " / " + std::to_string(BUFFER_CAP));
                
                for (auto const& [mod, t] : moduleLastTime) {
                    bool qEmpty = fileQueues[mod].empty();
                    std::string status = "Mod " + std::to_string(mod) + ": LastTime=" + std::to_string(t);
                    if (qEmpty) status += " [Queue Empty]";
                    else status += " [Queue OK]";
                    
                    if (t > currentMinTime) status += " [WAITING (Ahead)]";
                    else if (qEmpty) status += " [DONE]";
                    else status += " [STUCK?]";
                    
                    printLog(status);
                }
                
                // 強制脱出を試みる場合（デバッグ時のみ有効化）
                // forced = true ...
            }
        } else {
            // 何かしら動いていればカウンタをリセット
            stallCounter = 0;
        }

        // デッドロック回避 (Forced Read)
        if (stallCounter > 50000) { // 5万回空転したら強制的に読む
             printLog("[WARNING] Force reading to break deadlock...");
             for (auto& [m, q] : fileQueues) {
                if (!q.empty()) {
                    auto res = readEventsFromFile(q.front(), m, eventChunk, currentFileOffsets_[m]);
                    if (res.first > 0) moduleLastTime[m] = res.first;
                    if (res.second) { 
                         printLog("[INFO] Mod " + std::to_string(m) + ": Finished (Forced) " + q.front());
                         q.pop_front(); 
                         currentFileOffsets_[m] = 0; 
                    }
                    break; 
                }
            }
            stallCounter = 0; // リセット
        }
    }

END_ANALYSIS:
    printLog("Generating Hex Dump for debugging...");
    for (auto const& [modID, queue] : fileQueues) {
        if (!queue.empty()) printFileTail(queue.front(), currentFileOffsets_[modID]);
    }
    printLog("Analysis Loop Finished.");
}

// ----------------------------------------------------------------------------
// イベント処理 (画面更新 & ステータス管理)
// ----------------------------------------------------------------------------
bool DetectorAnalyzer::processChunk(const std::vector<Event>& sortedEvents) {
    if (sortedEvents.empty()) return true;

    // 定期的な画面更新 (1秒ごと)
    static auto lastUpdate = std::chrono::system_clock::now();
    auto now = std::chrono::system_clock::now();
    if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastUpdate).count() > 1000) {
        printSearchStatus();
        lastUpdate = now;
        dashboardShown_ = true;
    }
    if (!dashboardShown_) { printSearchStatus(); dashboardShown_ = true; }

    for (size_t i = 0; i < sortedEvents.size(); ++i) {
        const auto& e = sortedEvents[i];
        
        if (currentProcessingRunID_ != -1 && e.runID != currentProcessingRunID_) {
            printLog("Run transition detected. Resetting timing state.");
            currentRunStartTime_ = 0; 
            lastTimeoutCheckTime_ = 0;
            channelLastRefHitTime_.clear();
            moduleAliveTime_.clear(); // モジュール時計もクリア
            prevEventTimeForEff_ = 0; 
        }
        currentProcessingRunID_ = e.runID;

        // Run Start処理
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

        // Live Time計算 (モジュール独立時計に基づく判定)
        bool isSystemHealthy = true;
        for (int mid : activeModuleIDs_) {
            // モジュール時刻がデータ内にあるかチェック
            if (moduleAliveTime_.find(mid) == moduleAliveTime_.end()) {
                isSystemHealthy = false; break;
            }
        }
        if (isSystemHealthy && prevEventTimeForEff_ > 0) {
            unsigned long long dt_ns = e.eventTime_ns - prevEventTimeForEff_;
            // 異常な時間差(5秒以上)はスキップ
            if (dt_ns < 5000000000ULL) totalEffectiveTimeNs_ += dt_ns;
        }
        prevEventTimeForEff_ = e.eventTime_ns;

        // ゲイン解析トリガー
        if (lastGlobalAnalysisTime_ns_ == 0) lastGlobalAnalysisTime_ns_ = e.eventTime_ns;
        if (e.eventTime_ns >= lastGlobalAnalysisTime_ns_ + GAIN_ANALYSIS_WINDOW_NS) {
            analyzeGainShift(e.eventTime_ns);
            lastGlobalAnalysisTime_ns_ = e.eventTime_ns;
        }

        std::pair<int, int> kp = {e.module, e.sysCh};
        if (!ignoredMonitorChannels_.count(kp)) {
            channelLastRefHitTime_[e.module] = e.eventTime_ns;
        }
        foundChannels_[kp].insert(e.strip); // 生存ストリップ登録 (ラッチ)

        // ヒストグラムフィル
        if (e.module < MAX_MODULES && e.sysCh < MAX_SYS_CH) {
            if (ch_LUT[e.module][e.sysCh].pIndex >= 0) hGlobalPIndexMap->Fill(ch_LUT[e.module][e.sysCh].pIndex);
            if (hGainCheck_LUT[e.module][e.sysCh]) hGainCheck_LUT[e.module][e.sysCh]->Fill((double)e.tot);
            if (hRawToT_LUT[e.module][e.sysCh]) hRawToT_LUT[e.module][e.sysCh]->Fill((double)e.tot);
        }

        b_type = e.type; b_strip = e.strip; b_time = (long long)e.eventTime_ns; b_tot = e.tot;
        tree_->Fill();

        if (i > 0) {
            long long dt = (long long)e.eventTime_ns - (long long)sortedEvents[i-1].eventTime_ns;
            if (dt >= 0 && dt < 1000000) hDeltaT_Nearest->Fill((double)dt);
        }
    }
    return true;
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

void DetectorAnalyzer::determineIntegrationRange(TH1F* hist, int peakBin, int& outMinBin, int& outMaxBin) {
    int nBins = hist->GetNbinsX();
    const int SILENCE_BINS_5US = 250;
    
    outMinBin = 1;
    int zeroCount = 0;
    for (int i = peakBin; i > 1; --i) {
        if (hist->GetBinContent(i) <= 0) zeroCount++;
        else zeroCount = 0;
        if (zeroCount >= SILENCE_BINS_5US) {
            outMinBin = (i + zeroCount) - SILENCE_BINS_5US;
            if (outMinBin < 1) outMinBin = 1;
            break;
        }
    }

    outMaxBin = nBins;
    zeroCount = 0;
    for (int i = peakBin; i < nBins; ++i) {
        if (hist->GetBinContent(i) <= 0) zeroCount++;
        else zeroCount = 0;
        if (zeroCount >= SILENCE_BINS_5US) {
            outMaxBin = (i - zeroCount) + SILENCE_BINS_5US;
            if (outMaxBin > nBins) outMaxBin = nBins;
            break; 
        }
    }
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

int DetectorAnalyzer::findBestShift(TH1F* hTarget, TH1F* hTemplate, int binMin, int binMax, const std::string& /*debugLabel*/) {
    if (!hTarget || !hTemplate) return SHIFT_CALC_ERROR;
    if (hTarget->Integral(binMin, binMax) < 10) return SHIFT_CALC_ERROR; 
    
    double minResidual = 1e20;
    int bestShiftBins = SHIFT_CALC_ERROR;
    bool found = false;

    for (int s = -MATCHING_SEARCH_RANGE_BINS; s <= MATCHING_SEARCH_RANGE_BINS; ++s) {
        double res = calculateResidual(hTarget, hTemplate, s, binMin, binMax);
        if (res < 0) continue; 
        if (res < minResidual) { minResidual = res; bestShiftBins = s; found = true; }
    }
    if (!found) return SHIFT_CALC_ERROR;
    if (std::abs(bestShiftBins) == MATCHING_SEARCH_RANGE_BINS) return SHIFT_CALC_ERROR;
    return bestShiftBins;
}

bool DetectorAnalyzer::analyzeGainShift(unsigned long long currentTime_ns) { 
    if (!gainLogCsv_.is_open()) return true;

    if (!isCsvHeaderWritten_) {
        struct SortInfo { std::pair<int,int> k; int t; int s; std::string n; };
        std::vector<SortInfo> sorted;
        for (auto const& [keyPair, info] : detConfigMap_) sorted.push_back({keyPair, info.detTypeID, info.strip, info.cachedName});
        std::sort(sorted.begin(), sorted.end(), [](const SortInfo& a, const SortInfo& b){
            if(a.t != b.t) return a.t < b.t; return a.s < b.s;
        });
        gainLogCsv_ << "Time_ns"; rateLogCsv_ << "Time_ns";
        for(const auto& s : sorted) {
            std::string label = s.n + "-" + std::to_string(s.s);
            gainLogCsv_ << "," << label; rateLogCsv_ << "," << label;
        }
        gainLogCsv_ << "\n"; rateLogCsv_ << "\n";
        isCsvHeaderWritten_ = true;
    }

    gainLogCsv_ << currentTime_ns; rateLogCsv_ << currentTime_ns;
    
    std::vector<std::pair<int, int>> sortedKeys;
    for (auto const& [key, info] : detConfigMap_) sortedKeys.push_back(key);
    std::sort(sortedKeys.begin(), sortedKeys.end(), [&](std::pair<int,int> a, std::pair<int,int> b){
        const auto& iA = detConfigMap_[a]; const auto& iB = detConfigMap_[b];
        if(iA.detTypeID != iB.detTypeID) return iA.detTypeID < iB.detTypeID;
        return iA.strip < iB.strip; 
    });

    if (!isTemplateCaptured_) {
        for (const auto& key : sortedKeys) {
            int m = key.first; int c = key.second;
            cumulativeShiftMap_[key] = 0.0; 
            if (hGainCheck_LUT[m][c] && hGainCheck_LUT[m][c]->GetEntries() > 50) {
                int peakBin = findRightMostPeak(hGainCheck_LUT[m][c]);
                int rMin = 0, rMax = 0;
                determineIntegrationRange(hGainCheck_LUT[m][c], peakBin, rMin, rMax);
                rangeMinMap_[key] = rMin; rangeMaxMap_[key] = rMax;
                std::string tName = std::string(hGainCheck_LUT[m][c]->GetName()) + "_tmpl";
                if (templateHists_.count(key)) delete templateHists_[key];
                templateHists_[key] = (TH1F*)hGainCheck_LUT[m][c]->Clone(tName.c_str());
                templateHists_[key]->SetDirectory(0);
            }
            gainLogCsv_ << ",0";
            double rate = 0;
            if (hGainCheck_LUT[m][c]) {
                int rMin = rangeMinMap_[key]; int rMax = rangeMaxMap_[key];
                if (rMin > 0 && rMax > rMin) rate = hGainCheck_LUT[m][c]->Integral(rMin, rMax) / (double)(GAIN_ANALYSIS_WINDOW_NS / 1.0e9);
                hGainCheck_LUT[m][c]->Reset();
            }
            rateLogCsv_ << "," << rate;
        }
        isTemplateCaptured_ = true;
    } 
    else {
        for (const auto& key : sortedKeys) {
            int m = key.first; int c = key.second;
            std::string csvVal = "ND";
            double peakRate = 0.0;
            
            if (hGainCheck_LUT[m][c]) {
                TH1F* hist = hGainCheck_LUT[m][c];
                if (templateHists_.count(key) && hist->GetEntries() > 50) {
                    int rMin = rangeMinMap_[key]; int rMax = rangeMaxMap_[key];
                    std::string label = detConfigMap_[key].cachedName;
                    int deltaBins = findBestShift(hist, templateHists_[key], rMin, rMax, label);
                    
                    if (deltaBins != SHIFT_CALC_ERROR) {
                        rangeMinMap_[key] += deltaBins; rangeMaxMap_[key] += deltaBins;
                        double deltaNS = (double)deltaBins * BIN_WIDTH_NS;
                        cumulativeShiftMap_[key] += deltaNS;
                        csvVal = std::to_string(cumulativeShiftMap_[key]);

                        delete templateHists_[key];
                        std::string tName = std::string(hist->GetName()) + "_tmpl";
                        templateHists_[key] = (TH1F*)hist->Clone(tName.c_str());
                        templateHists_[key]->SetDirectory(0);
                        
                        peakRate = hist->Integral(rangeMinMap_[key], rangeMaxMap_[key]) / (double)(GAIN_ANALYSIS_WINDOW_NS / 1.0e9);
                    } else {
                        peakRate = hist->Integral(rMin, rMax) / (double)(GAIN_ANALYSIS_WINDOW_NS / 1.0e9);
                    }
                } else {
                    if (rangeMinMap_[key] > 0) peakRate = hist->Integral(rangeMinMap_[key], rangeMaxMap_[key]) / (double)(GAIN_ANALYSIS_WINDOW_NS / 1.0e9);
                }
                
                if (csvVal != "ND" && gainEvolutionGraphs_.count(key)) {
                    TGraph* g = gainEvolutionGraphs_[key];
                    g->SetPoint(g->GetN(), (double)currentTime_ns, cumulativeShiftMap_[key]);
                }
                if (rateEvolutionGraphs_.count(key)) {
                    TGraph* gr = rateEvolutionGraphs_[key];
                    gr->SetPoint(gr->GetN(), (double)currentTime_ns, peakRate);
                }
                hist->Reset();
            }
            gainLogCsv_ << "," << csvVal; rateLogCsv_ << "," << peakRate;
        }
    }
    gainLogCsv_ << std::endl; rateLogCsv_ << std::endl;
    return true; 
}

// ----------------------------------------------------------------------------
// ユーティリティ
// ----------------------------------------------------------------------------
void DetectorAnalyzer::printSearchStatus() {
    auto now = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - analysisStartTime_).count();
    
    double progress = 0.0;
    if (analysisDuration_ns_ != std::numeric_limits<unsigned long long>::max() && totalEffectiveTimeNs_ > 0) {
        double targetSec = (double)analysisDuration_ns_ / 1.0e9;
        double currentEffectiveSec = (double)totalEffectiveTimeNs_ / 1.0e9;
        double estimatedTotalBytes = (double)processedDataSize_ * (targetSec / currentEffectiveSec);
        if (estimatedTotalBytes > 0) progress = ((double)processedDataSize_ / estimatedTotalBytes) * 100.0;
    } else if (totalDataSize_ > 0) {
        progress = (double)processedDataSize_ / (double)totalDataSize_ * 100.0;
    }
    if (progress > 100.0) progress = 100.0;
    if (progress < 0.0) progress = 0.0;

    int barWidth = 40;
    int pos = (int)(barWidth * (progress / 100.0));

    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << " Progress: [";
    for (int i = 0; i < barWidth; ++i) std::cout << (i < pos ? "=" : (i == pos ? ">" : " "));
    std::cout << "] " << std::fixed << std::setprecision(1) << progress << " %" << std::endl;
    
    std::cout << " Live Time: " << (double)totalEffectiveTimeNs_ / 1.0e9 << " s";
    if (analysisDuration_ns_ != std::numeric_limits<unsigned long long>::max()) {
        std::cout << " / " << (int)(analysisDuration_ns_ / 1.0e9) << " s";
    }
    std::cout << " (Clock: " << (int)elapsed / 60 << "m " << (int)elapsed % 60 << "s)" << std::endl;

    // ★詳細なストリップ状態表示 (復活版: カラー + ラッチ)
    std::map<int, std::set<int>> activeStripsByDet;
    for (auto const& [key, strips] : foundChannels_) {
        int mod = key.first;
        int sys = key.second;
        if (ch_LUT[mod][sys].isValid) {
            int detID = ch_LUT[mod][sys].detTypeID;
            int str = ch_LUT[mod][sys].strip;
            activeStripsByDet[detID].insert(str);
        }
    }

    std::string detNames[] = {"X1", "Y1", "X2", "Y2"};
    std::cout << " Active Channels (Green = Detected):" << std::endl;

    for (int d = 0; d < 4; ++d) {
        std::cout << "  " << detNames[d] << ": ";
        const auto& activeSet = activeStripsByDet[d];
        for (int s = 1; s <= 32; ++s) { // 32chまで表示
            if (activeSet.count(s)) {
                std::cout << "\033[32m" << std::setw(2) << s << "\033[0m ";
            } else {
                std::cout << "\033[90m..\033[0m ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

std::string DetectorAnalyzer::getRunSignature(const std::string& fileName) {
    std::string fname = fileName;
    size_t lastSlash = fileName.find_last_of("/\\");
    if (lastSlash != std::string::npos) fname = fileName.substr(lastSlash + 1);
    static std::regex run_re(R"((PSD\d+_\d+))");
    std::smatch match;
    if (std::regex_search(fname, match, run_re)) return match[1].str();
    return "";
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
    for (const auto& info : allFiles) {
        fileQueues[info.modID].push_back(info.path);
    }
    for (auto const& [mod, q] : fileQueues) {
        printLog("Module " + std::to_string(mod) + ": " + std::to_string(q.size()) + " files registered (Sorted by Date->Seq).");
    }
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