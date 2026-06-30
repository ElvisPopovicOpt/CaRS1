#include <logger.hpp>
#include <iostream>
#include <filesystem>
#include <cassert>

namespace aco 
{

Logger::Logger(const std::string& logFile) 
    : logFilePath_(logFile)
{
    if (!logFile.empty()) 
    {
        // Kreiraj logs folder ako ne postoji
        std::filesystem::path logPath(logFile);
        if (logPath.has_parent_path()) 
        {
            std::filesystem::create_directories(logPath.parent_path());
        }
        
        // Otvori fajl za append (da se ne prepiše ako se pokrene više puta)
        logFile_.open(logFile, std::ios::app);
        if (!logFile_.is_open()) 
        {
            std::cerr << "Warning: Could not open log file: " << logFile << "\n";
        }
    }
}

Logger::~Logger() 
{
    if (logFile_.is_open()) 
    {
        logFile_.close();
    }
}

void Logger::log(const std::string& msg) 
{
    std::lock_guard<std::mutex> lock(mutex_);
    
    writeToConsole_(msg);
    if (logFile_.is_open()) 
    {
        writeToFile_(msg);
    }
}

void Logger::logNoNewline(const std::string& msg) 
{
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::cout << msg;
    std::cout.flush();
    
    if (logFile_.is_open()) 
    {
        logFile_ << msg;
        logFile_.flush();
    }
}

void Logger::flush() 
{
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::cout.flush();
    if (logFile_.is_open()) 
    {
        logFile_.flush();
    }
}

void Logger::writeToFile_(const std::string& msg) 
{
    if (logFile_.is_open()) 
    {
        logFile_ << msg << "\n";
        logFile_.flush();
    }
}

void Logger::writeToConsole_(const std::string& msg) 
{
    std::cout << msg << "\n";
    std::cout.flush();
}

} // namespace aco
