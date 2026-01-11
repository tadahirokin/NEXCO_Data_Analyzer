#include "DetectorAnalyzer.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <fstream> 
#include <stdexcept> 
#include <cstdio> 
#include <sys/stat.h> 
#include <cstring> 
#include <chrono>  
#include <ctime>   
#include <iomanip> 
#include <filesystem> // C++17の機能
#include <regex>

namespace fs = std::filesystem; // 名前空間の短縮

const int SYSTEM_CH_SIZE = 256; 
const int Y_CHANNEL_OFFSET = 64; 

void printLog(const std::string& msg) {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm* now_tm = std::localtime(&now_c);
    std::cout << "[" << std::put_time(now_tm, "%Y-%m-%d %H:%M:%S") << "] [INFO] " << msg << std::endl;
}

// Constructor
DetectorAnalyzer::DetectorAnalyzer(int timeWindow_ns, const std::string& outputFileName)
    : timeWindow_ns_(timeWindow_ns), 
      outputFile_(new TFile(outputFileName.c_str(), "RECREATE")),
      tree_(nullptr), hDeltaT_Nearest(nullptr), hDeltaT_n8(nullptr), hGlobalPIndexMap(nullptr),
      current_base_time_ns_(0),
      analysisDuration_ns_(std::numeric_limits<unsigned long long>::max()), 
      globalStartTimestamp_(0),
      next_gain_analysis_time_ns_(0), 
      isFirstGainAnalysis_(true),     
      isCsvHeaderWritten_(false)      
{
    if (!outputFile_ || outputFile_->IsZombie()) {
        throw std::runtime_error("[ERROR] Failed to create output ROOT file.");
    }
    setupTree();

    // フォルダ作成も filesystem を使うと楽ですが、念のため既存のまま
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

DetectorAnalyzer::~DetectorAnalyzer() {
    if (next_gain_analysis_time_ns_ > 0) {
        analyzeGainShift(next_gain_analysis_time_ns_);
    }

    if (outputFile_ && outputFile_->IsOpen()) {
        outputFile_->cd();
        if (tree_) tree_->Write();
        if (hDeltaT_Nearest) hDeltaT_Nearest->Write();
        if (hDeltaT_n8) hDeltaT_n8->Write();

        if (hGlobalPIndexMap) hGlobalPIndexMap->Write();

        struct HistSortEntry {
            int type;
            int strip;
            TH1F* hist;
        };
        
        std::vector<HistSortEntry> sortedRawHists;
        for (auto const& [key, hist] : hRawToTMap) {
            if (!hist) continue;
            int t = 999, s = key; 
            if (uniqueKeyToDetInfo_.count(key)) {
                t = uniqueKeyToDetInfo_[key].detTypeID;
                s = uniqueKeyToDetInfo_[key].detInCh;
            }
            sortedRawHists.push_back({t, s, hist});
        }

        std::sort(sortedRawHists.begin(), sortedRawHists.end(), 
            [](const HistSortEntry& a, const HistSortEntry& b) {
                if (a.type != b.type) return a.type < b.type;
                return a.strip < b.strip;
            }
        );

        for (const auto& entry : sortedRawHists) {
            entry.hist->Write();
        }

        std::vector<HistSortEntry> sortedGainHists;
        for (auto const& [key, hist] : gainCheckHists_) {
            if (!hist) continue;
            int t = 999, s = key;
            if (uniqueKeyToDetInfo_.count(key)) {
                t = uniqueKeyToDetInfo_[key].detTypeID;
                s = uniqueKeyToDetInfo_[key].detInCh;
            }
            sortedGainHists.push_back({t, s, hist});
        }

        std::sort(sortedGainHists.begin(), sortedGainHists.end(), 
            [](const HistSortEntry& a, const HistSortEntry& b) {
                if (a.type != b.type) return a.type < b.type;
                return a.strip < b.strip;
            }
        );

        for (const auto& entry : sortedGainHists) {
            entry.hist->Write();
        }

        outputFile_->Close();
        printLog("Analysis completed. Data saved (Sorted by Detector/Strip).");
    }
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "   Unmapped Channel Report (Threshold: >10 cps)" << std::endl;
    std::cout << "==========================================" << std::endl;

    double durationSec = 1.0;
    if (globalStartTimestamp_ > 0 && current_base_time_ns_ > globalStartTimestamp_) {
        durationSec = (double)(current_base_time_ns_ - globalStartTimestamp_) / 1.0e9;
    }
    if (durationSec < 1.0) durationSec = 1.0;

    bool foundSuspicious = false;

    for (auto const& [key, count] : unknownKeyCounts_) {
        double rate = (double)count / durationSec;
        if (rate > 10.0) {
            int mod = key / SYSTEM_CH_SIZE;
            int sys = key % SYSTEM_CH_SIZE;
            std::cout << "[WARN] HIGH RATE UNKNOWN SIGNAL: Mod " << mod 
                      << " Sys " << sys 
                      << " -> Total " << count << " hits (" 
                      << std::fixed << std::setprecision(1) << rate << " cps)" << std::endl;
            foundSuspicious = true;
        }
    }

    if (!foundSuspicious) {
        std::cout << "No high-rate unmapped signals found." << std::endl;
    } else {
        std::cout << "-> Check if these SysCh match your missing strips!" << std::endl;
    }
    std::cout << "==========================================\n" << std::endl;

    for (auto& pair : gainCheckHists_) { delete pair.second; }
    if (gainLogCsv_.is_open()) gainLogCsv_.close();
    if (rateLogCsv_.is_open()) rateLogCsv_.close();
}

void DetectorAnalyzer::setTimeLimitMinutes(double minutes) {
    if (minutes <= 0) {
        analysisDuration_ns_ = std::numeric_limits<unsigned long long>::max();
        printLog("Time limit set to: UNLIMITED");
    } else {
        analysisDuration_ns_ = (unsigned long long)(minutes * 60.0 * 1000000000.0);
        printLog("Time limit set to: " + std::to_string(minutes) + " minutes (" + std::to_string(analysisDuration_ns_) + " ns)");
    }
}

void DetectorAnalyzer::setupTree() {
    outputFile_->cd();
    tree_ = new TTree("Events", "Raw Detector Hits");
    
    tree_->Branch("type",  &b_type,  "type/I");   
    tree_->Branch("strip", &b_strip, "strip/I");  
    tree_->Branch("time",  &b_time,  "time/L");   
    tree_->Branch("tot",   &b_tot,   "tot/I");    

    hDeltaT_Nearest = new TH1F("DeltaT_Nearest", "Delta T (Nearest Neighbor);Delta T [ns];Counts", 10000, 0, 1000000);
    hDeltaT_n8 = new TH1F("DeltaT_n_plus_7", "Delta T (N to N+7);Delta T [ns];Counts", 2000, 0, 200000);

    hGlobalPIndexMap = new TH1F("GlobalHitMap_PIndex", "Global Hit Count vs P-Index;P-index;Total Counts", 600, 0, 600);
}

const std::map<int, std::vector<long long>>& DetectorAnalyzer::getRawTotData() const { return rawTotDataVec_; }

void DetectorAnalyzer::loadConversionFactorsFromCSV(const std::string& csvFileName, const std::string& detName, int moduleID) {
    std::ifstream file(csvFileName);
    if (!file.is_open()) throw std::runtime_error("[ERROR] Failed to open CSV: " + csvFileName);
    int typeID = -1;
    if (detName == "X1") typeID = ID_X1; else if (detName == "Y1") typeID = ID_Y1;
    else if (detName == "X2") typeID = ID_X2; else if (detName == "Y2") typeID = ID_Y2;

    std::string line;
    std::getline(file, line); 
    int loadedCount = 0;
    
    printLog("Loading CSV: " + csvFileName + " for " + detName + " (Mod " + std::to_string(moduleID) + ")");

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

            int localCh = 0;
            if (row[1].find('N') != std::string::npos || row[1].find('n') != std::string::npos) {
                continue; 
            } else {
                localCh = std::stoi(row[1]);
            }

            int pIndex = -1;
            if (row.size() >= 3 && !row[2].empty()) {
                try {
                    pIndex = std::stoi(row[2]);
                } catch (...) { pIndex = -1; }
            }

            int key = moduleID * SYSTEM_CH_SIZE + sysCh;
            std::string nameStr = (typeID < 2 ? "X" : "Y") + std::to_string(typeID % 2 + 1);
            uniqueKeyToDetInfo_[key] = { typeID, localCh, nameStr, pIndex };
            loadedCount++;

        } catch (...) { }
    }
    printLog(" -> Registered " + std::to_string(loadedCount) + " channels.");
}

// =========================================================================
//  ファイルの検索・ソート・分類を一括で行う関数 (std::filesystem 版)
// =========================================================================
void DetectorAnalyzer::loadAndSortFiles(const std::string& directory, std::map<int, std::deque<std::string>>& fileQueues) {
    // 1. ソート用のヘルパー構造体
    struct FileSortKey {
        long long dateVal; // YYYYMMDD
        int psdVal;        // PSDxxxxxx の番号
        std::string path;

        // ソート順定義: 日付(昇順) -> PSD番号(昇順) -> パス(昇順)
        bool operator<(const FileSortKey& other) const {
            if (dateVal != other.dateVal) return dateVal < other.dateVal;
            if (psdVal != other.psdVal) return psdVal < other.psdVal;
            return path < other.path;
        }
    };

    // 2. ファイル収集 (std::filesystem使用)
    std::vector<std::string> allFiles;
    
    if (fs::exists(directory)) {
        for (const auto& entry : fs::recursive_directory_iterator(directory)) {
            if (entry.is_regular_file() && entry.path().extension() == ".edb") {
                allFiles.push_back(entry.path().string());
            }
        }
    }

    if (allFiles.empty()) {
        std::cerr << "[ERROR] No .edb files found in " << directory << std::endl;
        return;
    }

    // 3. ソート処理
    std::sort(allFiles.begin(), allFiles.end(), [](const std::string& a, const std::string& b) {
        auto extractKey = [](const std::string& path) -> FileSortKey {
            FileSortKey key = {0, 0, path};
            try {
                // ファイルパスの中から "PSD" + "数字" + "_" + "数字(日付)" のパターンを探す
                size_t psdPos = path.rfind("PSD");
                if (psdPos != std::string::npos) {
                    std::string sub = path.substr(psdPos);
                    std::regex re(R"(PSD(\d+)_(\d{8}))");
                    std::smatch match;
                    if (std::regex_search(sub, match, re)) {
                        key.psdVal = std::stoi(match[1].str());
                        key.dateVal = std::stoll(match[2].str());
                    }
                }
            } catch(...) {}
            return key;
        };
        return extractKey(a) < extractKey(b);
    });

    // 4. モジュール振り分け
    int countMod0 = 0;
    int countMod1 = 0;
    for (const auto& f : allFiles) {
        if (f.find("_00_") != std::string::npos) {
            fileQueues[0].push_back(f);
            countMod0++;
        } else if (f.find("_01_") != std::string::npos) {
            fileQueues[1].push_back(f);
            countMod1++;
        } else {
            if (f.find("/PSD") != std::string::npos && f.find("_00_") != std::string::npos) {
                 fileQueues[0].push_back(f); countMod0++;
            } else if (f.find("/PSD") != std::string::npos && f.find("_01_") != std::string::npos) {
                 fileQueues[1].push_back(f); countMod1++;
            }
        }
    }

    printLog("Files loaded and sorted from: " + directory);
    std::cout << "  Mod 0 Files: " << countMod0 << std::endl;
    if (!fileQueues[0].empty()) std::cout << "   (First: " << fileQueues[0].front() << ")" << std::endl;
    std::cout << "  Mod 1 Files: " << countMod1 << std::endl;
    if (!fileQueues[1].empty()) std::cout << "   (First: " << fileQueues[1].front() << ")" << std::endl;
}

std::string DetectorAnalyzer::parseRunPrefix(const std::string& fullPath) {
    size_t lastSlash = fullPath.find_last_of("/\\");
    std::string fileName = (lastSlash == std::string::npos) ? fullPath : fullPath.substr(lastSlash + 1);
    size_t extPos = fileName.find_last_of(".");
    if (extPos != std::string::npos) fileName = fileName.substr(0, extPos);
    return fileName;
}

short DetectorAnalyzer::getRunID(const std::string& prefix) {
    if (runPrefixToId_.count(prefix)) return runPrefixToId_[prefix];
    short newId = (short)runIdToPrefix_.size();
    runIdToPrefix_.push_back(prefix);
    runPrefixToId_[prefix] = newId;
    return newId;
}

void DetectorAnalyzer::analyzeGainShift(unsigned long long currentTime_ns) {
    if (!gainLogCsv_.is_open() || !rateLogCsv_.is_open()) return;
    
    if (!isCsvHeaderWritten_) {
        gainLogCsv_ << "Time_ns"; rateLogCsv_ << "Time_ns";
        for (auto const& [key, info] : uniqueKeyToDetInfo_) {
            std::string label = info.cachedName + "-" + std::to_string(info.detInCh);
            gainLogCsv_ << "," << label;
            rateLogCsv_ << "," << label;
        }
        gainLogCsv_ << "\n"; rateLogCsv_ << "\n";
        isCsvHeaderWritten_ = true;
    }

    gainLogCsv_ << currentTime_ns; rateLogCsv_ << currentTime_ns;
    double windowSec = (double)GAIN_ANALYSIS_WINDOW_NS / 1.0e9;
    static TSpectrum *s = new TSpectrum(20); 

    for (auto const& [key, info] : uniqueKeyToDetInfo_) {
        double peakPos = 0.0; double peakRate = 0.0;
        
        if (gainCheckHists_.count(key)) {
            TH1F* hist = gainCheckHists_[key];
            bool isReferenceSet = refPeakPositions_.count(key);

            if (hist->GetEntries() > 0) {
                if (isFirstGainAnalysis_ || !isReferenceSet) {
                    if (hist->GetEntries() > 10) { 
                        Int_t nFound = s->Search(hist, 2.0, "goff", 0.1); 
                        if (nFound > 0) {
                            double maxToT = -1.0; int bestPeakIdx = -1;
                            for (int i = 0; i < nFound; ++i) {
                                double x = s->GetPositionX()[i];
                                if (x > maxToT) { maxToT = x; bestPeakIdx = i; }
                            }
                            if (bestPeakIdx != -1) {
                                peakPos = s->GetPositionX()[bestPeakIdx];
                                peakRate = s->GetPositionY()[bestPeakIdx] / windowSec;
                                refPeakPositions_[key] = peakPos; 
                            }
                        } else {
                            int maxBin = hist->GetMaximumBin();
                            if (hist->GetXaxis()->GetBinCenter(maxBin) > 100) {
                                peakPos = hist->GetXaxis()->GetBinCenter(maxBin);
                                peakRate = hist->GetBinContent(maxBin) / windowSec;
                                refPeakPositions_[key] = peakPos; 
                            }
                        }
                    }
                } else {
                    double refPos = refPeakPositions_[key];
                    int binMin = hist->GetXaxis()->FindBin(refPos * 0.8);
                    int binMax = hist->GetXaxis()->FindBin(refPos * 1.2);
                    if (binMin < 1) binMin = 1; if (binMax > hist->GetNbinsX()) binMax = hist->GetNbinsX();
                    hist->GetXaxis()->SetRange(binMin, binMax);
                    int maxBin = hist->GetMaximumBin();
                    if (hist->GetBinContent(maxBin) >= 1) {
                        peakPos = hist->GetXaxis()->GetBinCenter(maxBin);
                        peakRate = hist->GetBinContent(maxBin) / windowSec;
                    } else { peakPos = refPos; }
                    hist->GetXaxis()->SetRange(0, 0); 
                }
            }
            hist->Reset(); 
        }
        gainLogCsv_ << "," << peakPos; rateLogCsv_ << "," << peakRate;
    }
    gainLogCsv_ << std::endl; rateLogCsv_ << std::endl;
    gainLogCsv_ << std::flush;
    rateLogCsv_ << std::flush;

    if (isFirstGainAnalysis_) {
        isFirstGainAnalysis_ = false;
        printLog("Initial Gain Analysis completed.");
    }
}

void DetectorAnalyzer::processChunk(const std::vector<Event>& sortedEvents) {
    if (sortedEvents.empty()) return;
    
    size_t n = sortedEvents.size();
    for (size_t i = 0; i < n; ++i) {
        const auto& e = sortedEvents[i];

        if (next_gain_analysis_time_ns_ == 0) next_gain_analysis_time_ns_ = e.eventTime_ns + GAIN_ANALYSIS_WINDOW_NS;
        if (e.eventTime_ns >= next_gain_analysis_time_ns_) {
            analyzeGainShift(next_gain_analysis_time_ns_ - GAIN_ANALYSIS_WINDOW_NS);
            next_gain_analysis_time_ns_ += GAIN_ANALYSIS_WINDOW_NS;
            if (e.eventTime_ns >= next_gain_analysis_time_ns_) next_gain_analysis_time_ns_ = e.eventTime_ns + GAIN_ANALYSIS_WINDOW_NS;
        }

        int key = e.uniqueKey; 
        
        if (uniqueKeyToDetInfo_.count(key)) {
            int pIdx = uniqueKeyToDetInfo_[key].pIndex;
            if (pIdx >= 0) {
                hGlobalPIndexMap->Fill(pIdx);
            }
        }

        if (gainCheckHists_.find(key) == gainCheckHists_.end()) {
            std::string label = "Unknown";
            if (uniqueKeyToDetInfo_.count(key)) {
                const auto& info = uniqueKeyToDetInfo_[key];
                label = info.cachedName + "_Strip" + std::to_string(info.detInCh);
            }
            gainCheckHists_[key] = new TH1F(Form("GCheck_%s", label.c_str()), Form("Gain Check %s", label.c_str()), 10000, 0, 20000); 
            gainCheckHists_[key]->SetDirectory(0); 
        }
        gainCheckHists_[key]->Fill((double)e.tot);

        if (hRawToTMap.find(key) == hRawToTMap.end()) {
            outputFile_->cd();
            std::string hName, hTitle;
            if (uniqueKeyToDetInfo_.count(key)) {
                const auto& info = uniqueKeyToDetInfo_[key];
                hName = "RawToT_" + info.cachedName + "_Strip" + std::to_string(info.detInCh);
                hTitle = "Raw ToT: " + info.cachedName + " Strip " + std::to_string(info.detInCh);
            } else {
                hName = "RawToT_Unknown_Key" + std::to_string(key);
                hTitle = "Raw ToT: Unknown Key " + std::to_string(key);
            }
            hRawToTMap[key] = new TH1F(hName.c_str(), hTitle.c_str(), 20000, 0, 20000);
        }
        hRawToTMap[key]->Fill((double)e.tot);
        rawTotDataVec_[key].push_back(e.tot);

        b_type   = e.type;
        b_strip  = e.strip;
        b_time   = (long long)e.eventTime_ns; 
        b_tot    = e.tot;
        tree_->Fill();

        if (i > 0) {
            long long dt = (long long)e.eventTime_ns - (long long)sortedEvents[i-1].eventTime_ns;
            if (dt >= 0 && dt < 1000000) hDeltaT_Nearest->Fill((double)dt);
        }
        if (i >= 7) {
            long long dt8 = (long long)e.eventTime_ns - (long long)sortedEvents[i-7].eventTime_ns;
            if (dt8 >= 0 && dt8 < 200000) hDeltaT_n8->Fill((double)dt8);
        }
    }
}

std::pair<unsigned long long, bool> DetectorAnalyzer::readEventsFromFile(
    const std::string& fileName, 
    std::vector<Event>& rawEvents, 
    long long& offset
) {
    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs.is_open()) return {0, true}; 

    ifs.seekg(0, std::ios::end);
    long long fileSize = ifs.tellg();
    ifs.seekg(offset, std::ios::beg);

    std::string currentRunPrefix = parseRunPrefix(fileName);
    short currentRunID = getRunID(currentRunPrefix);
    
    const size_t BUFFER_SIZE = 64 * 1024 * 1024; 
    std::vector<char> buffer(BUFFER_SIZE);
    
    const size_t MAX_EVENTS_IN_THIS_CALL = 2000000; 
    size_t eventsReadThisCall = 0;

    size_t leftover = 0;
    unsigned long long lastEventTime = 0;
    bool timeLimitReached = false;
    
    static bool hasSeenTimeHeader = false; 

    while (ifs) {
        ifs.read(buffer.data() + leftover, BUFFER_SIZE - leftover);
        size_t readCount = ifs.gcount();
        if (readCount == 0 && leftover < 8) break; 

        size_t totalBytesInBuffer = leftover + readCount;
        size_t i = 0;
        unsigned char* buf = reinterpret_cast<unsigned char*>(buffer.data());

        while (i + 8 <= totalBytesInBuffer) {
            unsigned char header = buf[i];
            
            if (header == 0x69) { 
                unsigned char* p = &buf[i+1];
                unsigned int s = (static_cast<unsigned int>(p[0]) << 24) |
                                 (static_cast<unsigned int>(p[1]) << 16) |
                                 (static_cast<unsigned int>(p[2]) << 8)  |
                                 (static_cast<unsigned int>(p[3]));
                unsigned int ss = (static_cast<unsigned int>(p[4]) << 2) |
                                  ((static_cast<unsigned int>(p[5]) & 0xC0) >> 6);
                current_base_time_ns_ = (unsigned long long)s * 1000000000ULL + 
                                        (unsigned long long)ss * 1000000ULL;
                hasSeenTimeHeader = true; 
                i += 8; 
            } 
            else if (header == 0x6a) { 
                if (!hasSeenTimeHeader) { i++; continue; }
                unsigned char* p = &buf[i+1];
                unsigned int tof = (static_cast<unsigned int>(p[0]) << 16) |
                                   (static_cast<unsigned int>(p[1]) << 8)  |
                                   (static_cast<unsigned int>(p[2]));
                unsigned int pw  = (static_cast<unsigned int>(p[3]) << 12) |
                                   (static_cast<unsigned int>(p[4]) << 4)  |
                                   ((static_cast<unsigned int>(p[5]) & 0xF0) >> 4);
                unsigned int det = ((static_cast<unsigned int>(p[5]) & 0x0F) << 8) |
                                   (static_cast<unsigned int>(p[6]));
                int mod = (det >> 8) & 0xF;
                int sys = det & 0xFF;
                
                int key = mod * SYSTEM_CH_SIZE + sys;

                if (uniqueKeyToDetInfo_.count(key)) {
                    long long diff = static_cast<long long>(tof) - static_cast<long long>(pw);
                    unsigned long long t_ns = current_base_time_ns_ + diff;
                    lastEventTime = t_ns;
                    if (globalStartTimestamp_ == 0 || t_ns < globalStartTimestamp_) {
                        globalStartTimestamp_ = t_ns;
                    }
                    if (analysisDuration_ns_ != std::numeric_limits<unsigned long long>::max() &&
                        t_ns > globalStartTimestamp_ + analysisDuration_ns_) {
                        timeLimitReached = true;
                        i += 8; 
                        break; 
                    }
                    const auto& info = uniqueKeyToDetInfo_[key];
                    rawEvents.push_back({
                        info.detTypeID, info.detInCh, mod, t_ns, (int)pw, key
                    });
                    eventsReadThisCall++;
                } else {
                    unknownKeyCounts_[key]++;
                }
                i += 8; 
            } 
            else { i++; }
        }
        offset += i;
        if (timeLimitReached) return {lastEventTime, true}; 
        if (eventsReadThisCall >= MAX_EVENTS_IN_THIS_CALL) return {lastEventTime, false}; 
        leftover = totalBytesInBuffer - i;
        if (leftover > 0) std::memmove(buffer.data(), buffer.data() + i, leftover);
        if (offset >= fileSize) break;
    }
    return {lastEventTime, true}; 
}

void DetectorAnalyzer::processBinaryFiles(std::map<int, std::deque<std::string>>& fileQueues) {
    const size_t MAX_BUFFER_SIZE = 10000000; 
    std::vector<Event> eventChunk;
    eventChunk.reserve(MAX_BUFFER_SIZE); 
    std::map<int, unsigned long long> moduleLastTime;
    currentFileOffsets_.clear();
    for (auto const& [mod, queue] : fileQueues) {
        moduleLastTime[mod] = 0;
        currentFileOffsets_[mod] = 0;
    }
    auto compareFunc = [](const Event& a, const Event& b) { return a.eventTime_ns < b.eventTime_ns; };
    printLog("Starting Multi-Module Streaming Analysis...");

    while (true) {
        int targetMod = -1;
        unsigned long long minTime = std::numeric_limits<unsigned long long>::max();
        bool anyFileRemaining = false;
        for (auto& [mod, queue] : fileQueues) {
            if (queue.empty()) continue;
            if (moduleLastTime[mod] > 0 && globalStartTimestamp_ > 0 &&
                moduleLastTime[mod] > globalStartTimestamp_ + analysisDuration_ns_) {
                queue.clear();
                continue;
            }
            anyFileRemaining = true;
            if (moduleLastTime[mod] < minTime) {
                minTime = moduleLastTime[mod];
                targetMod = mod;
            }
        }
        if (!anyFileRemaining || targetMod == -1) break;
        std::string currentFile = fileQueues[targetMod].front();
        long long& offset = currentFileOffsets_[targetMod];
        auto result = readEventsFromFile(currentFile, eventChunk, offset);
        unsigned long long lastT = result.first;
        bool isFileDone = result.second;
        if (lastT > 0) moduleLastTime[targetMod] = lastT;
        if (isFileDone) {
            printLog("Finished file: " + currentFile);
            fileQueues[targetMod].pop_front();
            currentFileOffsets_[targetMod] = 0;
        }
        if (eventChunk.size() >= MAX_BUFFER_SIZE) {
            std::sort(eventChunk.begin(), eventChunk.end(), compareFunc);
            unsigned long long safeTime = std::numeric_limits<unsigned long long>::max();
            for (auto const& [mod, t] : moduleLastTime) {
                if (!fileQueues[mod].empty() && t < safeTime) safeTime = t;
            }
            if (safeTime == std::numeric_limits<unsigned long long>::max()) safeTime = eventChunk.back().eventTime_ns + 1; 
            size_t cutIndex = 0;
            for (size_t k = eventChunk.size(); k > 0; --k) {
                if (eventChunk[k-1].eventTime_ns <= safeTime) { 
                    cutIndex = k;
                    break;
                }
            }
            if (cutIndex > 0) {
                std::vector<Event> safeChunk(eventChunk.begin(), eventChunk.begin() + cutIndex);
                processChunk(safeChunk);
                std::vector<Event> remainingChunk;
                remainingChunk.reserve(MAX_BUFFER_SIZE);
                remainingChunk.insert(remainingChunk.end(), eventChunk.begin() + cutIndex, eventChunk.end());
                eventChunk = std::move(remainingChunk);
            }
        }
    }
    if (!eventChunk.empty()) {
        std::sort(eventChunk.begin(), eventChunk.end(), compareFunc);
        processChunk(eventChunk);
    }
}