#include "../include/CSVHandler.h"
#include <iostream>
#include <fstream>
#include <sstream>

// Constructor that accepts the file path
CSVHandler::CSVHandler(const std::string& filename) : filename_(filename) {}

// Function to read and parse the CSV file
std::vector<CSVRecord> CSVHandler::readCSV() {
    std::vector<CSVRecord> records;
    std::ifstream file(filename_);

    if (!file.is_open()) {
        std::cerr << "Could not open the file: " << filename_ << std::endl;
        return records;
    }

    std::string line;
    // Skip the header line
    std::getline(file, line);

    // Read each line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        CSVRecord record;

        try {
            // Parse sample number
            std::getline(ss, item, ',');
            record.sample_number = std::stoi(item);

            // Parse datetime
            std::getline(ss, item, ',');
            record.datetime = item;

            // Parse current value
            std::getline(ss, item, ',');
            record.curr_value = std::stod(item);

            // Parse voltage value
            std::getline(ss, item, ',');
            record.volt_value = std::stod(item);

            // Parse SoC value
            std::getline(ss, item, ',');
            record.soc_value = std::stod(item);

            // Add record to the list
            records.push_back(record);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing line: " << line << ". Exception: " << e.what() << "\n";
        }
    }

    file.close();
    return records;
}

// Function to save time series and SoC data to a CSV file (two-parameter version)
void CSVHandler::saveToCSV(const std::vector<double>& time_series, const std::vector<double>& soc) const {
    if (time_series.empty() || soc.empty()) {
        std::cerr << "Error: No data to write to file: " << filename_ << std::endl;
        return;
    }

    std::ofstream file(filename_);
    if (!file) {
        std::cerr << "Error opening file: " << filename_ << std::endl;
        return;
    }

    file << "Time,SoC\n";
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

// Function to save time series, SoC, current, and voltage data to a CSV file (four-parameter version)
void CSVHandler::saveToCSV(const std::vector<double>& time_series, const std::vector<double>& soc,
                           const std::vector<double>& current, const std::vector<double>& voltage) const {
    if (time_series.empty() || soc.empty() || current.empty() || voltage.empty()) {
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
        file << time_series[i] << "," << soc[i] << "," << current[i] << "," << voltage[i] << "\n";
        if (file.fail()) {
            std::cerr << "Error writing to file: " << filename_ << std::endl;
            file.close();
            return;
        }
    }

    file.close();
}
