#ifndef CSVREADER_H
#define CSVREADER_H

#include <string>
#include <vector>

// A structure to hold a single record from the CSV
struct CSVRecord {
    int sample_number;
    std::string datetime;
    double soc_value;
};

class CSVHandler {
public:
    // Constructor that accepts the file path
    CSVHandler(const std::string& filename);

    // Function to read and parse the CSV file
    std::vector<CSVRecord> readCSV();

    // Function to save time series and SoC data to a CSV file
    void saveToCSV(const std::vector<double>& time_series, const std::vector<double>& soc) const;

private:
    std::string filename_;
};

#endif // CSVREADER_H
