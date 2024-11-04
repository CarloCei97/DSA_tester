#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <cstdlib> // Include for system() function
//#include <matplot/matplot.h>
#include "include/CSVHandler.h"
#include "include/DownsamplingAlgorithms.h"
#include "include/WaveletTransform.h"
namespace fs = std::filesystem;
//namespace plt = matplot;

//###############################################################################################
//##################### handling time format conversion #########################################
//###############################################################################################

// Function to parse a datetime string and convert it to a time_point
std::chrono::system_clock::time_point parseDateTime(const std::string& datetime_str) {
    std::tm tm = {};
    std::istringstream ss(datetime_str);
    ss >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");
    auto tp = std::chrono::system_clock::from_time_t(std::mktime(&tm));

    // Handle fractional seconds (nanoseconds)
    size_t pos = datetime_str.find('.');
    if (pos != std::string::npos) {
        std::string fractional_seconds = datetime_str.substr(pos + 1);
        long nanoseconds = std::stol(fractional_seconds);

        // Convert time_point to a duration and add nanoseconds
        tp += std::chrono::duration_cast<std::chrono::system_clock::duration>(std::chrono::nanoseconds(nanoseconds));
    }

    return tp;
}

// Function to convert a time_point to a double (Unix timestamp)
double timePointToDouble(const std::chrono::system_clock::time_point& tp) {
    auto duration = tp.time_since_epoch();
    return std::chrono::duration_cast<std::chrono::seconds>(duration).count();
}

//###############################################################################################
//###############################################################################################
//###############################################################################################

int main() {
    //general parameters
    const int n_out = 2000;
    //douglas algorithm
    const double variance_threshold = 10;
    const double epsilon = 0.5;
    //wavelet transform
    std::string waveletType = "Haar";
    double threshold = 0.5;


    // Read data from CSV file
    std::string folder_path = "/Users/carlocei/Desktop/DS_A_Tester/real_data/";
    std::vector<std::string> csv_files;

    // List CSV files in the folder
    for (const auto& entry : fs::directory_iterator(folder_path)) {
        if (entry.path().extension() == ".csv") {
            csv_files.push_back(entry.path().filename().string());
        }
    }

    // Display the list of files for the user to choose from
    std::cout << "Available CSV files:\n";
    for (size_t i = 0; i < csv_files.size(); ++i) {
        std::cout << i + 1 << ": " << csv_files[i] << "\n";
    }

    // Get user's choice
    int choice;
    std::cout << "Enter the number of the file to select: ";
    std::cin >> choice;

    // Ensure valid choice
    if (choice < 1 || choice > csv_files.size()) {
        std::cerr << "Invalid choice.\n";
        return 1;
    }

    // Read data from the selected CSV file
    std::string selected_file = folder_path + csv_files[choice - 1];
    CSVHandler reader(selected_file);
    std::vector<CSVRecord> records = reader.readCSV();

// Vectors to store each column's data
    std::vector<double> time_series;
    std::vector<double> soc;
    std::vector<double> current;
    std::vector<double> voltage;

    for (const auto& record : records) {
        // Convert 'datetime' string to time_point, then to double
        auto time_point = parseDateTime(record.datetime); // Convert string to time_point
        double time_value = timePointToDouble(time_point); // Convert time_point to double

        // Push values into corresponding vectors
        time_series.push_back(time_value);
        soc.push_back(record.soc_value);
        current.push_back(record.curr_value);    // Assuming record has a field for current
        voltage.push_back(record.volt_value);    // Assuming record has a field for voltage
    }
    // Perform downsampling using DWT and keeping track of the execution time
    auto start_time = std::chrono::high_resolution_clock::now();
    auto [wavelet_sampled_time, wavelet_sampled_soc] = WaveletTransform::downsample(time_series, soc, waveletType, threshold);
    auto wavelet_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    // Perform downsampling using LTTB and keeping track of the execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [lttb_sampled_time, lttb_sampled_soc] = DownsamplingAlgorithms::largestTriangleThreeBuckets(time_series, soc, n_out);
    auto lttb_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();


    // Perform downsampling using Min-Max and keeping track of the execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [MinMax_sampled_time, MinMax_sampled_soc] = DownsamplingAlgorithms::minmaxDownsampling(time_series, soc, n_out);
    auto minmax_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

   // Perform downsampling using Piecewise Aggregate Approximation and keeping track of the execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [PAA_sampled_time, PAA_sampled_soc] = DownsamplingAlgorithms::piecewiseAggregateApproximation(time_series, soc, n_out);
    auto paa_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    // Perform downsampling using PLAA and keeping track of the execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [PLAA_sampled_time, PLAA_sampled_soc_withslope] = DownsamplingAlgorithms::piecewiseLinearAggregateApproximation(time_series, soc, n_out);
    auto PLAA_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
    // Vectors to hold the split values
    std::vector<double>  PLAA_sampled_soc;
    std::vector<double> PLAA_sampled_soc_slope;
    // Iterate through the original vector and split into two vectors
    for (const auto& pair : PLAA_sampled_soc_withslope) {
        PLAA_sampled_soc.push_back(pair.first);  // Add first element of pair to first_values
        PLAA_sampled_soc_slope.push_back(pair.second);  // Add second element of pair to second_values
    }

    // Perform downsampling using Random Sampling and keeping track of the execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [Random_sampled_time, Random_sampled_soc] = DownsamplingAlgorithms::randomSampling(time_series, soc, n_out);
    auto random_sampling_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    // Perform downsampling using Douglas-Peucker algorithm and keeping track of the execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [RDP_sampled_time, RDP_sampled_soc] = DownsamplingAlgorithms::RDPsampling(time_series, soc, epsilon);
    auto RDP_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();



    // Save the downsampled data to CSV files
    CSVHandler waveletHandler("../log/DWT_sampled.csv");
    CSVHandler lttbHandler("../log/lttb_sampled.csv");
    CSVHandler minmaxHandler("../log/minmax_sampled.csv");
    CSVHandler paaHandler("../log/paa_sampled.csv");
    CSVHandler adaptivePaaHandler("../log/PLAA_sampled.csv");
    CSVHandler randomHandler("../log/random_sampled.csv");
    CSVHandler hybridPaaMinMaxHandler("../log/hybrid_paa_minmax_sampled.csv");
    CSVHandler RDPHandler("../log/RDP_sampled.csv");
    waveletHandler.saveToCSV(wavelet_sampled_time, wavelet_sampled_soc);
    lttbHandler.saveToCSV(lttb_sampled_time, lttb_sampled_soc);
    minmaxHandler.saveToCSV(MinMax_sampled_time, MinMax_sampled_soc);
    paaHandler.saveToCSV(PAA_sampled_time, PAA_sampled_soc);
    adaptivePaaHandler.saveToCSV(PLAA_sampled_time, PLAA_sampled_soc);
    randomHandler.saveToCSV(Random_sampled_time, Random_sampled_soc);
    RDPHandler.saveToCSV(RDP_sampled_time, RDP_sampled_soc);

    // Plot the signals

//    // Original Signal
//    plt::figure();
//    plt::plot(time_series, soc);
//    plt::title("Original Signal");
//    plt::show();
//
//    // PAA downsampled
//    plt::figure();
//    plt::plot(time_series, Douglas_sampled);
//    plt::title("Douglas sampled Signal");
//    plt::show();

    // print the logfile
    std::ofstream logFile("/Users/carlocei/Desktop/DS_A_Tester/log/log.txt");
    logFile << "Discrete Wavelet Transform duration: " << wavelet_time << " seconds\n";
    logFile << "Largest Triangle Three Buckets duration: " << lttb_time << " seconds\n";
    logFile << "Min-Max Downsampling duration: " << minmax_time << " seconds\n";
    logFile << "Piecewise Aggregate Approximation duration: " << paa_time << " seconds\n";
    logFile << "Random Sampling duration: " << random_sampling_time << " seconds\n";
    logFile << "Adaptive PAA duration: " << PLAA_time << " seconds\n";
    logFile << "Douglas sampling duration: " << RDP_time << " seconds\n";
    logFile.close();

    // Launch the Python script to plot the graphs
    std::string pythonScript = "/Users/carlocei/Desktop/DS_A_Tester/plotgraphs.py";  // Full path to script.py
    std::string path = "/Users/carlocei/Desktop/DS_A_Tester/real_data/sample_data_SOC.csv";

    // Launch Python script with path as an argument
    std::string command = "python3 " + pythonScript + " " + path;
    int result = std::system(command.c_str());

    if (result != 0) {
        std::cerr << "Failed to execute Python script\n";
    }

    return 0;
}
