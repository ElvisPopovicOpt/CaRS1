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
    // Creates a logger that writes to console and optionally to a file.
    // If logFile is empty, only writes to console.
    explicit Logger(const std::string& logFile = "");

    ~Logger();

    // Writes a message to console and file (thread-safe)
    void log(const std::string& msg);

    // Writes a message without a trailing newline (thread-safe)
    void logNoNewline(const std::string& msg);

    // Flushes buffers (thread-safe)
    void flush();

    // True if the logger is usable (file is open, or no file was requested)
    bool isValid() const { return logFile_.is_open() || logFilePath_.empty(); }

private:
    std::string logFilePath_; // empty if not logging to a file
    std::ofstream logFile_;
    mutable std::mutex mutex_;
    
    void writeToFile_(const std::string& msg);
    void writeToConsole_(const std::string& msg);
};

} // namespace aco
