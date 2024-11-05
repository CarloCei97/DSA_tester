#ifndef CSVHANDLER_H
#define CSVHANDLER_H

#include <string>
#include <vector>
#include <chrono>

// Struct to hold individual CSV records
struct Record {
    int sampleNumber;
    std::string time;  // Original timestamp as string
    double timeValue;  // Converted timestamp as double (Unix time)
    double curr;
    double volt;
    double soc;
};

// Class to handle CSV file operations
class CSVHandler {
private:
    std::vector<Record> records;
    std::string filename_;

public:
    // Default constructor
    CSVHandler() = default;

    // Constructor that takes a filename
    CSVHandler(const std::string &filename) : filename_(filename) {}

    bool readCSV(const std::string &filePath);
    const std::vector<Record>& getRecords() const;

    // Function to save time series and SoC data to a CSV file
    void saveToCSV(const std::vector<double>& time_series, const std::vector<double>& soc) const;

    // Static methods for time conversion
    static std::chrono::system_clock::time_point parseDateTime(const std::string& datetime_str);
    static double timePointToDouble(const std::chrono::system_clock::time_point& tp);
};

#endif // CSVHANDLER_H
