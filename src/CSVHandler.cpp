#include "../include/CSVHandler.h"
#include <fstream>
#include <sstream>
#include <iostream>

bool CSVHandler::readCSV(const std::string &filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return false;
    }

    std::string line;
    bool isFirstLine = true;

    while (std::getline(file, line)) {
        if (isFirstLine) {
            isFirstLine = false;
            continue;
        }

        std::stringstream ss(line);
        std::string cell;
        Record record;

        std::getline(ss, cell, ',');  // Sample number
        record.sampleNumber = std::stoi(cell);

        std::getline(ss, record.time, ',');  // Timestamp
        std::getline(ss, cell, ',');  // Current
        record.curr = std::stod(cell);
        std::getline(ss, cell, ',');  // Voltage
        record.volt = std::stod(cell);
        std::getline(ss, cell, ',');  // SOC
        record.soc = std::stod(cell);

        records.push_back(record);
    }

    file.close();
    return true;
}

const std::vector<Record>& CSVHandler::getRecords() const {
    return records;
}

void CSVHandler::saveToCSV(const std::vector<double>& time_series, const std::vector<double>& soc ) const {
    if (time_series.empty() || soc.empty() ) {
        std::cerr << "Error: No data to write to file: " << filename_ << std::endl;
        return;
    }

    std::ofstream file(filename_);
    if (!file) {
        std::cerr << "Error opening file: " << filename_ << std::endl;
        return;
    }

    file << "Time,SoC,Current,Voltage\n";
    for (size_t i = 0; i < time_series.size(); ++i) {
        file << time_series[i] << "," << soc[i] << "\n";
        if (file.fail()) {
            std::cerr << "Error writing to file: " << filename_ << std::endl;
            file.close();
            return;
        }
    }

    file.close();
}
