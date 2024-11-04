#ifndef CSVREADER_H
#define CSVREADER_H

#include <string>
#include <vector>

// A structure to hold a single record from the CSV
struct CSVRecord {
    int sample_number;
    std::string datetime;
    double soc_value;
    double curr_value;
    double volt_value;
};

class CSVHandler {
public:
    explicit CSVHandler(const std::string& filename);

    // Function to read CSV records
    std::vector<CSVRecord> readCSV();

    // Overloaded function to save only time_series and soc
    void saveToCSV(const std::vector<double>& time_series, const std::vector<double>& soc) const;

    // Overloaded function to save time_series, soc, current, and voltage
    void saveToCSV(const std::vector<double>& time_series, const std::vector<double>& soc,
                   const std::vector<double>& current, const std::vector<double>& voltage) const;

private:
    std::string filename_;
};


#endif // CSVREADER_H
