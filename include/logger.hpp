#pragma once
#include <string>
#include <fstream>
#include <mutex>
#include <memory>

namespace aco 
{

class Logger 
{
public:
    // Kreira logger koji zapisuje u konzolu i u fajl
    // Ako je logFile prazan, zapisuje samo u konzolu
    explicit Logger(const std::string& logFile = "");
    
    ~Logger();
    
    // Zapisuje poruku i u konzolu i u fajl (thread-safe)
    void log(const std::string& msg);
    
    // Zapisuje poruku bez newline (thread-safe)
    void logNoNewline(const std::string& msg);
    
    // Flush buffer (thread-safe)
    void flush();
    
    // Provjeri je li logger validan (fajl je otvoren ili je logFilePath_ prazan)
    bool isValid() const { return logFile_.is_open() || logFilePath_.empty(); }
    
private:
    std::string logFilePath_; // prazan string ako se ne koristi fajl
    std::ofstream logFile_;
    mutable std::mutex mutex_;
    
    void writeToFile_(const std::string& msg);
    void writeToConsole_(const std::string& msg);
};

} // namespace aco
