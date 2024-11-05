#ifndef CSVHANDLER_H
#define CSVHANDLER_H

#include <string>
#include <vector>

// Struct to hold individual CSV records
struct Record {
    int sampleNumber;
    std::string time;
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
};

#endif // CSVHANDLER_H
