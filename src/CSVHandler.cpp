#include "../include/CSVHandler.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

// Function to read data from a CSV file
bool CSVHandler::readCSV(const std::string &filePath) {
    // Open the file for reading
    std::ifstream file(filePath);
    if (!file.is_open()) {
        // Error message if the file could not be opened
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return false;
    }

    std::string line;
    bool isFirstLine = true;  // Flag to skip the header line

    // Read each line of the file
    while (std::getline(file, line)) {
        // Skip the header
        if (isFirstLine) {
            isFirstLine = false;
            continue;
        }

        std::stringstream ss(line);
        std::string cell;
        Record record;

        // Parse each field from the CSV line

        std::getline(ss, cell, ',');  // Read sample number
        record.sampleNumber = std::stoi(cell);  // Convert to integer

        std::getline(ss, record.time, ',');  // Read timestamp as a string
        record.timeValue = timePointToDouble(parseDateTime(record.time));  // Convert timestamp to Unix time

        std::getline(ss, cell, ',');  // Read current value
        record.curr = std::stod(cell);  // Convert to double

        std::getline(ss, cell, ',');  // Read voltage value
        record.volt = std::stod(cell);  // Convert to double

        std::getline(ss, cell, ',');  // Read SOC value
        record.soc = std::stod(cell);  // Convert to double

        // Add the parsed record to the records vector
        records.push_back(record);
    }

    // Close the file after reading
    file.close();
    return true;
}

// Function to return the list of parsed records
const std::vector<Record>& CSVHandler::getRecords() const {
    return records;
}

// Function to save time series and SoC data to a CSV file
void CSVHandler::saveToCSV(const std::vector<double>& time_series, const std::vector<double>& soc) const {
    // Check if there is data to write
    if (time_series.empty() || soc.empty()) {
        std::cerr << "Error: No data to write to file: " << filename_ << std::endl;
        return;
    }

    // Open the file for writing
    std::ofstream file(filename_);
    if (!file) {
        std::cerr << "Error opening file: " << filename_ << std::endl;
        return;
    }

    // Write header
    file << "Time,SoC\n";

    // Write each data point to the file
    for (size_t i = 0; i < time_series.size(); ++i) {
        file << time_series[i] << "," << soc[i] << "\n";

        // Check if writing to file was successful
        if (file.fail()) {
            std::cerr << "Error writing to file: " << filename_ << std::endl;
            file.close();
            return;
        }
    }

    // Close the file after writing
    file.close();
}


//#######################################################################################################
//##################### handling time format conversion #################################################
//#######################################################################################################

// Static method to parse a datetime string and convert it to a time_point
std::chrono::system_clock::time_point CSVHandler::parseDateTime(const std::string& datetime_str) {
    std::tm tm = {};  // Initialize a time structure
    std::istringstream ss(datetime_str);
    ss >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");  // Parse date and time

    // Convert tm to a time_point
    auto tp = std::chrono::system_clock::from_time_t(std::mktime(&tm));

    // Handle fractional seconds (nanoseconds) if present
    size_t pos = datetime_str.find('.');
    if (pos != std::string::npos) {
        // Extract fractional seconds as nanoseconds
        std::string fractional_seconds = datetime_str.substr(pos + 1);
        long nanoseconds = std::stol(fractional_seconds);

        // Add nanoseconds to the time_point
        tp += std::chrono::duration_cast<std::chrono::system_clock::duration>(std::chrono::nanoseconds(nanoseconds));
    }

    return tp;
}

// Static method to convert a time_point to a double (Unix timestamp)
double CSVHandler::timePointToDouble(const std::chrono::system_clock::time_point& tp) {
    // Convert time_point to Unix time in seconds
    auto duration = tp.time_since_epoch();
    return std::chrono::duration_cast<std::chrono::seconds>(duration).count();
}

//#######################################################################################################
//#######################################################################################################
//#######################################################################################################
