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

bool DetectorAnalyzer::processChunk(const std::vector<Event>& sortedEvents) {
    if (enableMuonAnalysis_) {
        return processChunk_Muon(sortedEvents);
    } else {
        return processChunk_Fast(sortedEvents);
    }
}

// ----------------------------------------------------------------------------
// ★爆速モード: ゲイン解析のみ (バッファリング・ミュオン抽出なし)
// ----------------------------------------------------------------------------
bool DetectorAnalyzer::processChunk_Fast(const std::vector<Event>& sortedEvents) {
    if (sortedEvents.empty()) return true;
    static bool isReported[MAX_MODULES][MAX_SYS_CH] = {false};

    for (size_t i = 0; i < sortedEvents.size(); ++i) {
        const auto& e = sortedEvents[i];
        
        // 1. 基本解析 (ラン初期化)
        if (currentRunStartTime_ == 0) {
            if (e.tot >= 1000) {
                currentRunStartTime_ = e.eventTime_ns;
                lastTimeoutCheckTime_ = e.eventTime_ns; 
                prevEventTimeForEff_ = e.eventTime_ns;
                lastGlobalAnalysisTime_ns_ = e.eventTime_ns; // ここで初期化しておく
                for(int mid : activeModuleIDs_) moduleAliveTime_[mid] = e.eventTime_ns;
            } else continue;
        }
        currentEventTime_ns_ = e.eventTime_ns;
        moduleAliveTime_[e.module] = e.eventTime_ns;

        bool isSystemHealthy = true;
        for (int mid : activeModuleIDs_) {
            if (moduleAliveTime_.find(mid) == moduleAliveTime_.end()) { isSystemHealthy = false; break; }
        }
        
        // processBinaryFilesで上書きされるため実質無効だが、互換性維持のため残す
        if (isSystemHealthy && prevEventTimeForEff_ > 0) {
            unsigned long long dt_ns = e.eventTime_ns - prevEventTimeForEff_;
            if (dt_ns < 5000000000ULL) totalEffectiveTimeNs_ += dt_ns;
        }
        prevEventTimeForEff_ = e.eventTime_ns;

        // 2. データの充填
        if (e.module < MAX_MODULES && e.sysCh < MAX_SYS_CH) {
            if (!isReported[e.module][e.sysCh]) {
                foundChannels_[{e.module, e.sysCh}].insert(e.strip);
                isReported[e.module][e.sysCh] = true; 
            }
            // 解析用ヒストグラムに今のイベントを入れる
            if (hGainCheck_LUT[e.module][e.sysCh]) hGainCheck_LUT[e.module][e.sysCh]->Fill((double)e.tot);
            
            // 生ToTヒストグラムにも入れる
            if (hRawToT_LUT[e.module][e.sysCh]) hRawToT_LUT[e.module][e.sysCh]->Fill((double)e.tot);
            
            int pIdx = ch_LUT[e.module][e.sysCh].pIndex;
            if (hGlobalPIndexMap && pIdx >= 0) hGlobalPIndexMap->Fill((double)pIdx);
        }

        // ★削除済: イベント時刻依存の analyzeGainShift 呼び出しを削除 (SafeTimeベースへ移行)

        // 4. Tree Fill
        b_type = static_cast<Char_t>(e.type);
        b_strip = static_cast<Char_t>(e.strip);
        b_tot = static_cast<Int_t>(e.tot);
        b_time = static_cast<Long64_t>(e.eventTime_ns - globalStartTimestamp_);
        if (tree_) tree_->Fill();

        // 5. DeltaT
        if (i > 0) {
            long long dt1 = (long long)e.eventTime_ns - (long long)sortedEvents[i-1].eventTime_ns;
            if (dt1 >= 0 && dt1 < 1000000 && hDeltaT_Nearest) hDeltaT_Nearest->Fill((double)dt1);
        }
        if (i + 7 < sortedEvents.size()) {
            long long dt8 = (long long)sortedEvents[i+7].eventTime_ns - (long long)e.eventTime_ns;
            if (dt8 >= 0 && dt8 < 1000000 && hDeltaT_n8) hDeltaT_n8->Fill((double)dt8);
        }
    }
    return true;
}

// ----------------------------------------------------------------------------
// ★Slowモード: ゲイン解析 + ミュオン抽出
// ----------------------------------------------------------------------------
bool DetectorAnalyzer::processChunk_Muon(const std::vector<Event>& sortedEvents) {
    if (sortedEvents.empty()) return true;
    static bool isReported[MAX_MODULES][MAX_SYS_CH] = {false};

    for (size_t i = 0; i < sortedEvents.size(); ++i) {
        const auto& e = sortedEvents[i];
        
        // 1. 基本解析 (ラン初期化)
        if (currentRunStartTime_ == 0) {
            if (e.tot >= 1000) {
                currentRunStartTime_ = e.eventTime_ns;
                lastTimeoutCheckTime_ = e.eventTime_ns; 
                prevEventTimeForEff_ = e.eventTime_ns;
                lastGlobalAnalysisTime_ns_ = e.eventTime_ns; // 初期化
                for(int mid : activeModuleIDs_) moduleAliveTime_[mid] = e.eventTime_ns;
            } else continue;
        }
        currentEventTime_ns_ = e.eventTime_ns;
        moduleAliveTime_[e.module] = e.eventTime_ns;

        bool isSystemHealthy = true;
        for (int mid : activeModuleIDs_) {
            if (moduleAliveTime_.find(mid) == moduleAliveTime_.end()) { isSystemHealthy = false; break; }
        }
        if (isSystemHealthy && prevEventTimeForEff_ > 0) {
            unsigned long long dt_ns = e.eventTime_ns - prevEventTimeForEff_;
            if (dt_ns < 5000000000ULL) totalEffectiveTimeNs_ += dt_ns;
        }
        prevEventTimeForEff_ = e.eventTime_ns;

        // 2. データの充填
        if (e.module < MAX_MODULES && e.sysCh < MAX_SYS_CH) {
            if (!isReported[e.module][e.sysCh]) {
                foundChannels_[{e.module, e.sysCh}].insert(e.strip);
                isReported[e.module][e.sysCh] = true; 
            }
            // 解析用ヒストグラムに今のイベントを入れる
            if (hGainCheck_LUT[e.module][e.sysCh]) hGainCheck_LUT[e.module][e.sysCh]->Fill((double)e.tot);
            
            // 生ToTヒストグラムにも入れる
            if (hRawToT_LUT[e.module][e.sysCh]) hRawToT_LUT[e.module][e.sysCh]->Fill((double)e.tot);

            int pIdx = ch_LUT[e.module][e.sysCh].pIndex;
            if (hGlobalPIndexMap && pIdx >= 0) hGlobalPIndexMap->Fill((double)pIdx);
        }

        // ★削除済: イベント時刻依存の analyzeGainShift 呼び出しを削除 (SafeTimeベースへ移行)

        // 4. Tree Fill
        b_type = static_cast<Char_t>(e.type);
        b_strip = static_cast<Char_t>(e.strip);
        b_tot = static_cast<Int_t>(e.tot);
        b_time = static_cast<Long64_t>(e.eventTime_ns - globalStartTimestamp_);
        if (tree_) tree_->Fill();

        // 5. DeltaT
        if (i > 0) {
            long long dt1 = (long long)e.eventTime_ns - (long long)sortedEvents[i-1].eventTime_ns;
            if (dt1 >= 0 && dt1 < 1000000 && hDeltaT_Nearest) hDeltaT_Nearest->Fill((double)dt1);
        }
        if (i + 7 < sortedEvents.size()) {
            long long dt8 = (long long)sortedEvents[i+7].eventTime_ns - (long long)e.eventTime_ns;
            if (dt8 >= 0 && dt8 < 1000000 && hDeltaT_n8) hDeltaT_n8->Fill((double)dt8);
        }

        // Slowモード固有: バッファリング
        analysisBuffer_.push_back(e);
    }

    // Slowモード固有: ミュオン抽出実行
    processEventExtraction();

    return true;
}

void DetectorAnalyzer::processEventExtraction() {
    if (analysisBuffer_.empty()) return;

    while (analysisBuffer_.size() >= 8) { 
        struct Candidate {
            size_t startIdx;
            size_t length; 
            double width;
        };

        Candidate bestCand = {0, 0, 1.0e20}; 
        bool foundAny = false;
        
        size_t scanLimit = 8;
        if (analysisBuffer_.size() < scanLimit + 8) scanLimit = analysisBuffer_.size() - 8;

        for (size_t i = 0; i < scanLimit; ++i) {
            double w = 0;
            bool isA = checkConditionA(analysisBuffer_, i, w);
            bool isB = false;
            
            if (!isA) isB = checkConditionB(analysisBuffer_, i, w);

            if (isA || isB) {
                if (w < bestCand.width) {
                    bestCand = {i, (isA ? 8UL : 7UL), w};
                    foundAny = true;
                }
            }
        }

        if (foundAny) {
            hEventWidth_All->Fill(bestCand.width);
            if (bestCand.length == 8) hEventWidth_CondA->Fill(bestCand.width);
            else hEventWidth_CondB->Fill(bestCand.width);

            for (size_t k = 0; k < bestCand.startIdx; ++k) {
                analysisBuffer_.pop_front();
            }

            b_i_runID = currentProcessingRunID_;
            b_i_eventID = globalEventID_Ideal_++;
            b_i_nHits = (int)bestCand.length;
            b_i_width = bestCand.width;
            
            b_i_type->clear(); b_i_strip->clear(); b_i_tot->clear(); b_i_time->clear();
            
            double tMin = analysisBuffer_[0].eventTime_ns; 

            for (size_t k = 0; k < bestCand.length; ++k) {
                const auto& e = analysisBuffer_[k];
                b_i_type->push_back(ch_LUT[e.module][e.sysCh].detTypeID);
                b_i_strip->push_back(ch_LUT[e.module][e.sysCh].strip);
                b_i_tot->push_back(e.tot);
                b_i_time->push_back((double)(e.eventTime_ns - tMin));
            }
            if (idealTree_) idealTree_->Fill();

            for (size_t k = 0; k < bestCand.length; ++k) {
                analysisBuffer_.pop_front();
            }

        } else {
            analysisBuffer_.pop_front();
        }
    }
}

// --- 条件判定ヘルパー関数 ---
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

    for(int type=0; type<4; ++type) {
        if (!hasAdjacentPair(sMap[type])) return false;
    }
    return true;
}

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

    if (!hasAdjacentPair(sMap[0])) return false;
    if (!hasAdjacentPair(sMap[2])) return false;
    if (!hasAdjacentPair(sMap[3])) return false;
    if (sMap[1].size() != 1) return false;
    return isAdjacentToOffline(sMap[1][0], 1);
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

bool DetectorAnalyzer::syncDataStream(std::ifstream& ifs, long long& skippedBytes) {
    skippedBytes = 0;
    const size_t SCAN_BUFFER_SIZE = 1024 * 1024; 
    std::vector<char> buffer(SCAN_BUFFER_SIZE);
    
    long long globalOffset = ifs.tellg(); 
    long long initialPos = globalOffset;
    
    while (ifs) {
        ifs.read(buffer.data(), SCAN_BUFFER_SIZE);
        size_t bytesRead = ifs.gcount();
        if (bytesRead == 0) break;

        for (size_t i = 0; i < bytesRead; ++i) {
            unsigned char h = static_cast<unsigned char>(buffer[i]);
            if (h == 0x69 || h == 0x6a) {
                bool syncOK = true;
                if (i + 80 <= bytesRead) {
                    for (int k = 1; k < 10; ++k) {
                        unsigned char nextH = static_cast<unsigned char>(buffer[i + (k * 8)]);
                        if (nextH != 0x69 && nextH != 0x6a) {
                            syncOK = false;
                            break;
                        }
                    }
                } else {
                    long long candidatePos = globalOffset + i;
                    long long nextChunkPos = globalOffset + bytesRead; 
                    for (int k = 1; k < 10; ++k) {
                         ifs.seekg(candidatePos + (k * 8), std::ios::beg);
                         char nextH;
                         if (!ifs.read(&nextH, 1)) { syncOK = false; break; }
                         if ((unsigned char)nextH != 0x69 && (unsigned char)nextH != 0x6a) {
                             syncOK = false; break;
                         }
                    }
                    if (!syncOK) {
                        ifs.clear();
                        ifs.seekg(nextChunkPos, std::ios::beg);
                    }
                }

                if (syncOK) {
                    long long finalPos = globalOffset + i;
                    ifs.seekg(finalPos, std::ios::beg);
                    skippedBytes = finalPos - initialPos;
                    return true;
                }
            }
        }
        globalOffset += bytesRead;
        long long currentSkipped = globalOffset - initialPos;
        if (currentSkipped > 0 && currentSkipped % (10 * 1024 * 1024) == 0) {
             printLog("[SYNC] Scanning... Skipped: " + std::to_string(currentSkipped / 1024 / 1024) + " MB");
        }
    }
    skippedBytes = globalOffset - initialPos;
    return false;
}

void DetectorAnalyzer::setupTree() {
    outputFile_->cd();
    tree_ = new TTree("Events", "Raw Detector Hits");
    tree_->Branch("type",  &b_type,  "type/B");
    tree_->Branch("strip", &b_strip, "strip/B");
    tree_->Branch("tot",   &b_tot,   "tot/I");
    tree_->Branch("time",  &b_time,  "time/L");

    hGlobalPIndexMap = new TH1F("hGlobalPIndexMap", "Global Hit Count vs P-Index;p-index;Counts", 600, 0, 600);
    hDeltaT_Nearest = new TH1F("DeltaT_Nearest", "Delta T (Nearest Neighbor);Delta T [ns];Counts", 100000, 0, 100000);
    hDeltaT_n8      = new TH1F("DeltaT_n_plus_7", "Delta T (N to N+7);Delta T [ns];Counts", 100000, 0, 100000);

    if (enableMuonAnalysis_) {
        hEventWidth_All   = new TH1F("hEventWidth_All",   "Event Width (All Clusters);Width [ns];Counts", 200, 0, 1000);
        hEventWidth_CondA = new TH1F("hEventWidth_CondA", "Event Width (Cond A: 8-hit);Width [ns];Counts", 200, 0, 1000);
        hEventWidth_CondB = new TH1F("hEventWidth_CondB", "Event Width (Cond B: 7-hit);Width [ns];Counts", 200, 0, 1000);

        TDirectory* currentDir = gDirectory;
        idealFile_ = new TFile("ideal_events.root", "RECREATE");
        idealTree_ = new TTree("IdealEvents", "Muon Events (Shortest Width Selected)");

        b_i_type  = new std::vector<int>();
        b_i_strip = new std::vector<int>();
        b_i_tot   = new std::vector<int>();
        b_i_time  = new std::vector<double>();

        idealTree_->Branch("runID",   &b_i_runID,   "runID/I");
        idealTree_->Branch("eventID", &b_i_eventID, "eventID/I");
        idealTree_->Branch("nHits",   &b_i_nHits,   "nHits/I");
        idealTree_->Branch("width",   &b_i_width,   "width/D");
        idealTree_->Branch("type",  &b_i_type);
        idealTree_->Branch("strip", &b_i_strip);
        idealTree_->Branch("tot",   &b_i_tot);
        idealTree_->Branch("time",  &b_i_time);

        globalEventID_Ideal_ = 0;
        currentDir->cd();
    } else {
        hEventWidth_All = nullptr; hEventWidth_CondA = nullptr; hEventWidth_CondB = nullptr;
        idealFile_ = nullptr; idealTree_ = nullptr;
        b_i_type = nullptr; b_i_strip = nullptr; b_i_tot = nullptr; b_i_time = nullptr;
    }
}

void DetectorAnalyzer::generatePDFReport() {
    printLog("Generating PDF Report...");
    std::string pdfName = "AnalysisReport.pdf";
    TCanvas* c1 = new TCanvas("c1", "Report", 1200, 800);
    c1->Print((pdfName + "(").c_str()); 

    c1->SetLogy(1);
    if (hDeltaT_Nearest) { hDeltaT_Nearest->Draw(); c1->Print(pdfName.c_str()); }
    if (hDeltaT_n8)      { hDeltaT_n8->Draw();      c1->Print(pdfName.c_str()); }
    if (hEventWidth_All)   { hEventWidth_All->Draw();   c1->Print(pdfName.c_str()); }
    if (hEventWidth_CondA) { hEventWidth_CondA->Draw(); c1->Print(pdfName.c_str()); }
    if (hEventWidth_CondB) { hEventWidth_CondB->Draw(); c1->Print(pdfName.c_str()); }

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
        h->SetTitle(Form("Raw ToT (Log): %s Strip%d (Mod%d Ch%d)", detConfigMap_[key].cachedName.c_str(), detConfigMap_[key].strip, entry.mod, entry.ch));
        h->Draw();

        int rMin = rangeMinMap_[key];
        int rMax = rangeMaxMap_[key];
        
        if (rMin > 0 && rMax > rMin) {
            double x1 = 0; double x2 = 0;
            if (gainCheckHists_.find(key) != gainCheckHists_.end()) {
                TH1F* refHist = gainCheckHists_[key];
                if (refHist) {
                    x1 = refHist->GetXaxis()->GetBinLowEdge(rMin);
                    x2 = refHist->GetXaxis()->GetBinUpEdge(rMax);
                }
            } else {
                x1 = (double)rMin * 20.0;
                x2 = (double)rMax * 20.0;
            }
            double yMin = (c1->GetLogy()) ? 0.5 : 0;
            TBox* box = new TBox(x1, yMin, x2, h->GetMaximum());
            box->SetFillStyle(0); box->SetLineColor(kRed); box->SetLineWidth(2);
            box->Draw("same");
        }
        c1->Print(pdfName.c_str());
    }

    c1->SetLogy(0);
    c1->Print((pdfName + ")").c_str()); 
    delete c1;
}

DetectorAnalyzer::~DetectorAnalyzer() {
    if (currentEventTime_ns_ > lastGlobalAnalysisTime_ns_) analyzeGainShift(currentEventTime_ns_);
    calculateEffectiveTime();
    generatePDFReport();

    if (outputFile_ && outputFile_->IsOpen()) {
        outputFile_->cd();
        if (tree_) tree_->Write();
        if (hGlobalPIndexMap) hGlobalPIndexMap->Write();
        if (hEventWidth_All) hEventWidth_All->Write();
        if (hEventWidth_CondA) hEventWidth_CondA->Write();
        if (hEventWidth_CondB) hEventWidth_CondB->Write();

        auto dirHist = outputFile_->mkdir("Histograms"); 
        if (dirHist) {
            dirHist->cd();
            for (auto const& [key, h] : hRawToTMap) if (h) h->Write();
        }

        outputFile_->cd(); 
        auto dirTemp = outputFile_->mkdir("Templates"); 
        if (dirTemp) {
            dirTemp->cd();
            for (auto const& [key, h] : templateHists_) if (h) h->Write();
        }
        
        outputFile_->cd();
        auto dG = outputFile_->mkdir("GainHistory");
        if (dG) {
            dG->cd();
            for (auto const& [key, g] : gainEvolutionGraphs_) {
                if (g && g->GetN() > 0) g->Write();
            }
        }

        outputFile_->cd();
        auto dR = outputFile_->mkdir("RateHistory");
        if (dR) {
            dR->cd();
            for (auto const& [key, g] : rateEvolutionGraphs_) {
                if (g && g->GetN() > 0) g->Write();
            }
        }
        outputFile_->Close();
        delete outputFile_;
        outputFile_ = nullptr;
    }

    if (idealFile_ && idealFile_->IsOpen()) {
        idealFile_->cd();
        if (idealTree_) idealTree_->Write();
        idealFile_->Close();
        delete idealFile_;
        idealFile_ = nullptr;
        printLog("Ideal events saved to: ideal_events.root");
    }

    if (gainLogCsv_.is_open()) gainLogCsv_.close();
    if (rateLogCsv_.is_open()) rateLogCsv_.close();
}

// ----------------------------------------------------------------------------
// ファイル読み込み (T0アンカー + 異常値は論理同期解除)
// ----------------------------------------------------------------------------

std::pair<unsigned long long, bool> DetectorAnalyzer::readEventsFromFile(const std::string& fileName, int modID, std::vector<Event>& rawEvents, long long& offset) {
    if (offset == 0) {
        printLog("[OPENING] Mod " + std::to_string(modID) + ": " + fileName);
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
    bool hasFoundValidT0 = hasSeenTimeHeaderMap_[modID]; 
    
    // 新規ファイル時はT0探索モードへリセット
    if (offset == 0) hasFoundValidT0 = false;

    const size_t BUFFER_SIZE = 64 * 1024 * 1024; 
    std::vector<char> buffer(BUFFER_SIZE);
    unsigned long long lastT = currentBaseTime; 

    std::vector<size_t> t0_indices;           
    std::vector<size_t> t0_event_snapshots;   

    ifs.read(buffer.data(), BUFFER_SIZE);
    size_t readCount = ifs.gcount();
    
    if (readCount > 0) {
        processedDataSize_ += readCount; 
        size_t i = 0;
        unsigned char* buf = reinterpret_cast<unsigned char*>(buffer.data());

        while (i + 8 <= readCount) {
            unsigned char h = buf[i];
            bool isDataAligned = false;

            // ---------------------------------------------------------
            // 1. 物理的なデータ整列チェック (isDataAligned の確定)
            // ---------------------------------------------------------
            if (i + 80 <= readCount) {
                if (h == 0x69 || h == 0x6a) { 
                    isDataAligned = true; 
                    for (int k = 1; k < 10; ++k) {
                        unsigned char nh = buf[i + (k * 8)];
                        if (nh != 0x69 && nh != 0x6a) { isDataAligned = false; break; }
                    }
                }
            } else {
                // 末端処理: データらしければAlignedとみなす
                if (h == 0x69 || h == 0x6a) isDataAligned = true;
            }

            // ---------------------------------------------------------
            // 2. モード別処理
            // ---------------------------------------------------------
            if (!isDataAligned) {
                // 物理同期が取れていない -> 即座に1バイト探索へ
                // (T0が見つかっていようがいまいが、物理的にズレたらSync Lost)
                hasFoundValidT0 = false; 
                i++;
                continue;
            }

            // ここから先は isDataAligned == true (8バイト刻み処理) が確定
            
            if (!hasFoundValidT0) {
                // --- T0探索モード ---
                // データ区切りは正しい。あとは有効なT0(重複以降の時刻)を探す。
                if (h == 0x69) {
                    unsigned char* p = &buf[i+1];
                    unsigned long long s = (static_cast<unsigned long long>(p[0]) << 24) | (p[1] << 16) | (p[2] << 8) | p[3];
                    unsigned long long ss = (static_cast<unsigned long long>(p[4]) << 2) | ((p[5] & 0xC0) >> 6);
                    unsigned long long t = ((p[5] & 0x3F) << 2) | ((p[6] & 0xC0) >> 6);
                    unsigned long long foundTime = s * 1000000000ULL + ss * 1000000ULL + t;

                    if (foundTime > baseTimeMap_[modID]) {
                        long long diff = (long long)(foundTime - currentBaseTime);
                        if (currentBaseTime == 0 || std::abs(diff) < 1000000000LL) {
                            if (currentBaseTime != 0) {
                                deadTimeRanges_.push_back({currentBaseTime, foundTime});
                            }
                            currentBaseTime = foundTime;
                            hasFoundValidT0 = true; 

                            // T0確定時の諸処理
                            lastT = currentBaseTime;
                            moduleAliveTime_[modID] = currentBaseTime;
                            if (firstT0Time_ == 0) firstT0Time_ = s;
                            t0_indices.push_back(i);
                            t0_event_snapshots.push_back(rawEvents.size());
                        }
                    }
                }
                // T0が見つかるまでは読み飛ばし(8バイト)
                i += 8;
                
            } else {
                // --- 通常解析モード ---
                // T0確定済み。データをパースする。
                if (h == 0x69) { 
                    unsigned char* p = &buf[i+1];
                    unsigned long long s = (static_cast<unsigned long long>(p[0]) << 24) | (p[1] << 16) | (p[2] << 8) | p[3];
                    unsigned long long ss = (static_cast<unsigned long long>(p[4]) << 2) | ((p[5] & 0xC0) >> 6);
                    unsigned long long t = ((p[5] & 0x3F) << 2) | ((p[6] & 0xC0) >> 6);
                    unsigned long long newTime = s * 1000000000ULL + ss * 1000000ULL + t;

                    long long diff = (long long)(newTime - currentBaseTime);
                    
                    if (std::abs(diff) > 1000000000LL) {
                        hasFoundValidT0 = false; // Sync Lost
                    } else {
                        currentBaseTime = newTime;
                        lastT = newTime;
                        moduleAliveTime_[modID] = newTime;
                    }
                    t0_indices.push_back(i);
                    t0_event_snapshots.push_back(rawEvents.size());

                } else if (h == 0x6a) { 
                    unsigned char* p = &buf[i+1];
                    unsigned int tof = (p[0] << 16) | (p[1] << 8) | p[2];
                    unsigned int pw  = (p[3] << 12) | (p[4] << 4) | ((p[5] & 0xF0) >> 4);
                    int det = ((p[5] & 0x0F) << 8) | p[6];
                    int mod = (det >> 8) & 0xF;
                    int sys = det & 0xFF;
                    
                    if (mod < MAX_MODULES && sys < MAX_SYS_CH && ch_LUT[mod][sys].isValid) {
                        if (tof >= pw) { 
                            unsigned long long eventT = currentBaseTime + (unsigned long long)(tof - pw);
                            lastT = eventT;
                            moduleAliveTime_[modID] = eventT;
                            rawEvents.push_back({ ch_LUT[mod][sys].detTypeID, ch_LUT[mod][sys].strip, mod, eventT, (int)pw, sys });
                        }
                    }
                }
                i += 8;
            }
        }
        
        bool willBeEOF = (offset + readCount >= (unsigned long long)fileSize);
        if (!willBeEOF && t0_indices.size() >= 2) {
            size_t lastT0Pos = t0_indices.back();
            size_t snapshotSize = t0_event_snapshots.back();
            offset += lastT0Pos;
            if (rawEvents.size() > snapshotSize) {
                rawEvents.resize(snapshotSize);
            }
        } else {
            offset += readCount;
        }
        
    } else { offset = fileSize; }

    baseTimeMap_[modID] = currentBaseTime;
    hasSeenTimeHeaderMap_[modID] = hasFoundValidT0;
    return {lastT, (offset >= fileSize)}; 
}



bool DetectorAnalyzer::analyzeGainShift(unsigned long long currentTime_ns) {
    double timeMin = currentTime_ns / 60.0e9;

    if (!isCsvHeaderWritten_) {
        std::string detNames[] = {"X1", "Y1", "X2", "Y2"};
        gainLogCsv_ << "Time_ns";
        rateLogCsv_ << "Time_ns";
        for (int type = 0; type < 4; ++type) {
            for (int s = 1; s <= 32; ++s) {
                std::string header = "," + detNames[type] + "_" + std::to_string(s);
                gainLogCsv_ << header;
                rateLogCsv_ << header;
            }
        }
        gainLogCsv_ << "\n";
        rateLogCsv_ << "\n";
        isCsvHeaderWritten_ = true;
    }

    std::vector<double> currentStepPeaks(128, 0.0);
    std::vector<double> currentStepRates(128, 0.0);

    for (auto const& [key, hist] : gainCheckHists_) {
        auto& cfg = detConfigMap_[key];
        if (!hist) continue;

        if (templateHists_.find(key) == templateHists_.end()) {
            int peakBin = findRightMostPeak(hist);
            int rMin = 0, rMax = 0;
            
            if (peakBin > 0) {
                determineIntegrationRange(hist, peakBin, rMin, rMax);
            }

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
        else if (rangeMinMap_[key] > 0) {
            TH1F* tHist = templateHists_[key];
            int deltaBin = findBestShift(hist, tHist, rangeMinMap_[key], rangeMaxMap_[key], cfg.cachedName);
            
            if (deltaBin != SHIFT_CALC_ERROR) {
                rangeMinMap_[key] += deltaBin;
                rangeMaxMap_[key] += deltaBin;
                cumulativeShiftNsMap_[key] += (deltaBin * BIN_WIDTH_NS);

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
                currentStepRates[idx] = counts / (GAIN_ANALYSIS_WINDOW_NS / 1.0e9) * 60.0;
            }
        }
        hist->Reset();
    }

    gainLogCsv_ << currentTime_ns;
    rateLogCsv_ << currentTime_ns;
    for (int i = 0; i < 128; ++i) {
        gainLogCsv_ << "," << currentStepPeaks[i];
        rateLogCsv_ << "," << currentStepRates[i];
    }
    gainLogCsv_ << "\n";
    rateLogCsv_ << "\n";

    return true;
}

// ----------------------------------------------------------------------------
// ★修正: メインループ
// (SafeTime基準解析トリガー + 累積デッドタイム減算方式 + 5秒ルール)
// ----------------------------------------------------------------------------
void DetectorAnalyzer::processBinaryFiles(std::map<int, std::deque<std::string>>& fileQueues) {
    const size_t PER_MOD_CAP = 1000000; 
    std::map<int, std::deque<Event>> moduleBuffers;
    std::map<int, unsigned long long> moduleLastTime;
    auto lastUIDraw = std::chrono::system_clock::now();

    for (auto const& [mod, q] : fileQueues) {
        moduleLastTime[mod] = 0; currentFileOffsets_[mod] = 0;
        hasSeenTimeHeaderMap_[mod] = false; baseTimeMap_[mod] = 0;
    }
    deadTimeRanges_.clear();
    isAnalysisStarted_ = false; globalStartTimestamp_ = 0;
    
    std::string currentRunSignature = ""; 
    unsigned long long lastSafeTime = 0;
    unsigned long long accumulatedDeadTimeNs = 0;

    printSearchStatus();

    while (true) {
        if (currentRunSignature.empty()) {
            for (auto const& [m, q] : fileQueues) {
                if (!q.empty()) {
                    currentRunSignature = getRunSignature(q.front());
                    printLog("=== STARTING RUN: " + currentRunSignature + " ===");
                    currentRunStartTime_ = 0; lastGlobalAnalysisTime_ns_ = 0;
                    isAnalysisStarted_ = false; deadTimeRanges_.clear(); 
                    accumulatedDeadTimeNs = 0;
                    for(auto& [mod, t] : moduleLastTime) t = 0;
                    break; 
                }
            }
        }
        if (currentRunSignature.empty()) break; 

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
                lastGlobalAnalysisTime_ns_ = maxBase + ANALYSIS_START_SKIP_NS; 
                isAnalysisStarted_ = true;
                printLog("Sync Established. Gain Analysis scheduled to start after 3 mins skip.");
            }
        }

        if (isAnalysisStarted_) {
            unsigned long long safeTime = getSafeTime(moduleLastTime);
            if (safeTime > lastSafeTime) {
                unsigned long long diff = safeTime - lastSafeTime;
                
                if (diff >= 5000000000ULL) {
                    accumulatedDeadTimeNs += diff;
                }

                if (safeTime >= globalStartTimestamp_ + accumulatedDeadTimeNs) {
                     totalEffectiveTimeNs_ = (safeTime - globalStartTimestamp_) - accumulatedDeadTimeNs;
                }

                if (totalEffectiveTimeNs_ >= analysisDuration_ns_) {
                    printLog("[LIMIT] Time limit reached. Stopping analysis loop.");
                    break; 
                }

                lastSafeTime = safeTime; currentEventTime_ns_ = safeTime;
                
                if (safeTime >= lastGlobalAnalysisTime_ns_ + GAIN_ANALYSIS_WINDOW_NS) {
                    analyzeGainShift(safeTime);
                    lastGlobalAnalysisTime_ns_ = safeTime;
                }
                
                while(!deadTimeRanges_.empty() && deadTimeRanges_.front().second < safeTime) {
                    deadTimeRanges_.erase(deadTimeRanges_.begin());
                }
            }
            
            std::vector<Event> mergeChunk;
            for (auto& [m, buf] : moduleBuffers) {
                while (!buf.empty() && buf.front().eventTime_ns <= safeTime) {
                    Event ev = std::move(buf.front()); buf.pop_front();
                    
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

        bool runOver = true;
        for (auto const& [m, q] : fileQueues) if (!q.empty() && getRunSignature(q.front()) == currentRunSignature) runOver = false;
        for (auto const& [m, b] : moduleBuffers) if (!b.empty()) runOver = false;
        if (runOver) {
            if(currentEventTime_ns_ > lastGlobalAnalysisTime_ns_) analyzeGainShift(currentEventTime_ns_);
            currentRunSignature = "";
        }
    }
}

void DetectorAnalyzer::printSearchStatus() {
    static size_t lastReportedCount = 0;

    if (!isAnalysisStarted_) {
        std::cout << "\r\033[K" << "[STATUS] Syncing..." << std::flush;
        return;
    }

    auto now = std::chrono::system_clock::now();
    double elapsedClock = std::chrono::duration_cast<std::chrono::seconds>(now - analysisStartTime_).count();
    
    double liveTimeMin = (double)totalEffectiveTimeNs_ / 1.0e9 / 60.0;
    double progress = (analysisDuration_ns_ > 0) ? (double)totalEffectiveTimeNs_ / (double)analysisDuration_ns_ * 100.0 : 0.0;
    
    std::cout << "\r\033[K" 
              << "Progress: [" << std::fixed << std::setprecision(1) << progress << "%] "
              << "Data: " << liveTimeMin << " min | " 
              << "Real: " << (int)elapsedClock / 60 << "m" << (int)elapsedClock % 60 << "s | "
              << "Live: " << liveTimeMin << " min " 
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
    int binMin = hist->GetXaxis()->FindBin(10000.0); 
    int binMax = hist->GetXaxis()->FindBin(100000.0); 
    if (binMax > hist->GetNbinsX()) binMax = hist->GetNbinsX();
    if (binMin < 1) binMin = 1;

    for (int i = binMax; i >= binMin; --i) {
        if (hist->GetBinContent(i) > 0) {
            return i; 
        }
    }
    return 0; 
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

void DetectorAnalyzer::determineIntegrationRange(TH1F* hist, int peakBin, int& rMin, int& rMax) {
    if (!hist || peakBin <= 0) {
        rMin = 0; rMax = 0;
        return;
    }

    double binWidthNs = hist->GetXaxis()->GetBinWidth(1); 
    if (binWidthNs <= 0.0) binWidthNs = 20.0; 
    const int BINS_2US = (int)(2000.0 / binWidthNs); 
    const double MIN_WIDTH_NS = 4000.0; 

    rMax = std::min(hist->GetNbinsX(), peakBin + BINS_2US);
    int foundGapStart = 0;
    
    for (int i = peakBin - 1; i >= 1; --i) {
        bool isGap = true;
        for (int k = 0; k < BINS_2US; ++k) {
            if ((i - k) < 1) {
                isGap = true; 
                foundGapStart = 1;
                break;
            }
            if (hist->GetBinContent(i - k) > 0) {
                isGap = false;
                break;
            }
        }
        if (isGap) {
            foundGapStart = std::max(1, i - (BINS_2US - 1));
            break;
        }
    }

    if (foundGapStart == 0) {
        rMin = 0; rMax = 0;
        return;
    }
    rMin = foundGapStart;

    double currentWidthNs = (rMax - rMin) * binWidthNs;
    if (currentWidthNs <= MIN_WIDTH_NS) {
        rMin = 0; rMax = 0;
        return;
    }
}

int DetectorAnalyzer::findBestShift(TH1F* hTarget, TH1F* hTemplate, int binMin, int binMax, const std::string&) {
    double minRes = 1e20; int bestS = SHIFT_CALC_ERROR; bool found = false;
    const int SEARCH_RANGE = 25; 

    for (int s = -SEARCH_RANGE; s <= SEARCH_RANGE; ++s) { 
        double res = calculateResidual(hTarget, hTemplate, s, binMin, binMax);
        if (res >= 0 && res < minRes) { minRes = res; bestS = s; found = true; }
    }
    if (!found || std::abs(bestS) == SEARCH_RANGE) return SHIFT_CALC_ERROR;
    return bestS;
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

    for (const auto& info : allFiles) {
        fileQueues[info.modID].push_back(info.path);
    }

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