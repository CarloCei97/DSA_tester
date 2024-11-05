#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <cstdlib> //system()
#include "include/CSVHandler.h"
#include "include/DownsamplingAlgorithms.h"
#include "include/WaveletTransform.h"
namespace fs = std::filesystem;


int main() {
//########################################## parameters #################################################
    //general parameters
    const int n_out = 2000;
    //douglas algorithm
    const double variance_threshold = 10;
    const double epsilon = 0.5;
    //wavelet transform
    std::string waveletType = "Haar";
    double threshold = 0.5;
//##########################################################################################################

//########################################## import data from file #########################################
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
    CSVHandler csvHandler;
    csvHandler.readCSV(selected_file);
    const std::vector<Record>& records = csvHandler.getRecords();


// Vectors to store each column's data
    std::vector<double> time_series;
    std::vector<double> soc;
    std::vector<double> current;
    std::vector<double> voltage;

    for (const auto& record : records) {
        // Push values into corresponding vectors
        time_series.push_back(record.timeValue);
        soc.push_back(record.soc);
        current.push_back(record.curr);    // Assuming record has a field for current
        voltage.push_back(record.volt);    // Assuming record has a field for voltage
    }
//#######################################################################################################

//################################# perform the downsampling algorithms #################################

// Perform downsampling using LTTB for each variable separately and keep track of execution time
    auto start_time = std::chrono::high_resolution_clock::now();
    auto [lttb_sampled_time_soc, lttb_sampled_soc] = DownsamplingAlgorithms::largestTriangleThreeBuckets(time_series, soc, n_out);
    auto lttb_time_soc = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [lttb_sampled_time_volt, lttb_sampled_volt] = DownsamplingAlgorithms::largestTriangleThreeBuckets(time_series, voltage, n_out);
    auto lttb_time_volt = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [lttb_sampled_time_curr, lttb_sampled_curr] = DownsamplingAlgorithms::largestTriangleThreeBuckets(time_series, current, n_out);
    auto lttb_time_curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

// Perform downsampling using Min-Max for each variable separately and keep track of execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [MinMax_sampled_time_soc, MinMax_sampled_soc] = DownsamplingAlgorithms::minmaxDownsampling(time_series, soc, n_out);
    auto minmax_time_soc = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [MinMax_sampled_time_volt, MinMax_sampled_volt] = DownsamplingAlgorithms::minmaxDownsampling(time_series, voltage, n_out);
    auto minmax_time_volt = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [MinMax_sampled_time_curr, MinMax_sampled_curr] = DownsamplingAlgorithms::minmaxDownsampling(time_series, current, n_out);
    auto minmax_time_curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

// Perform downsampling using Piecewise Aggregate Approximation (PAA) for each variable separately and keep track of execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [PAA_sampled_time_soc, PAA_sampled_soc] = DownsamplingAlgorithms::piecewiseAggregateApproximation(time_series, soc, n_out);
    auto paa_time_soc = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [PAA_sampled_time_volt, PAA_sampled_volt] = DownsamplingAlgorithms::piecewiseAggregateApproximation(time_series, voltage, n_out);
    auto paa_time_volt = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [PAA_sampled_time_curr, PAA_sampled_curr] = DownsamplingAlgorithms::piecewiseAggregateApproximation(time_series, current, n_out);
    auto paa_time_curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    // Perform downsampling using Random Sampling for each variable separately and keep track of execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [Random_sampled_time_soc, Random_sampled_soc] = DownsamplingAlgorithms::randomSampling(time_series, soc, n_out);
    auto random_sampling_time_soc = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [Random_sampled_time_volt, Random_sampled_volt] = DownsamplingAlgorithms::randomSampling(time_series, voltage, n_out);
    auto random_sampling_time_volt = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [Random_sampled_time_curr, Random_sampled_curr] = DownsamplingAlgorithms::randomSampling(time_series, current, n_out);
    auto random_sampling_time_curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

// Perform downsampling using Douglas-Peucker (RDP) for each variable separately and keep track of execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [RDP_sampled_time_soc, RDP_sampled_soc] = DownsamplingAlgorithms::RDPsampling(time_series, soc, epsilon);
    auto RDP_time_soc = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [RDP_sampled_time_volt, RDP_sampled_volt] = DownsamplingAlgorithms::RDPsampling(time_series, voltage, epsilon);
    auto RDP_time_volt = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [RDP_sampled_time_curr, RDP_sampled_curr] = DownsamplingAlgorithms::RDPsampling(time_series, current, epsilon);
    auto RDP_time_curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    /*
// Perform downsampling using DWT for each variable separately and keep track of execution time
    auto start_time = std::chrono::high_resolution_clock::now();
    auto [wavelet_sampled_time_soc, wavelet_sampled_soc] = WaveletTransform::downsample(time_series, soc, waveletType, threshold);
    auto wavelet_time_soc = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    auto start_time = std::chrono::high_resolution_clock::now();
    auto [wavelet_sampled_time_volt, wavelet_sampled_volt] = WaveletTransform::downsample(time_series, voltage, waveletType, threshold);
    auto wavelet_time_volt = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto [wavelet_sampled_time_curr, wavelet_sampled_curr] = WaveletTransform::downsample(time_series, current, waveletType, threshold);
    auto wavelet_time_curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
*/
//#######################################################################################################

//############################ save the results of the downsampling #####################################

// LTTB downsampling files
    CSVHandler lttbHandler_soc("../log/lttb_sampled_soc.csv");
    CSVHandler lttbHandler_curr("../log/lttb_sampled_curr.csv");
    CSVHandler lttbHandler_volt("../log/lttb_sampled_volt.csv");

// MinMax downsampling files
    CSVHandler minmaxHandler_soc("../log/minmax_sampled_soc.csv");
    CSVHandler minmaxHandler_curr("../log/minmax_sampled_curr.csv");
    CSVHandler minmaxHandler_volt("../log/minmax_sampled_volt.csv");

// PAA downsampling files
    CSVHandler paaHandler_soc("../log/paa_sampled_soc.csv");
    CSVHandler paaHandler_curr("../log/paa_sampled_curr.csv");
    CSVHandler paaHandler_volt("../log/paa_sampled_volt.csv");

// Adaptive PAA downsampling files
    CSVHandler adaptivePaaHandler_soc("../log/PLAA_sampled_soc.csv");
    CSVHandler adaptivePaaHandler_curr("../log/PLAA_sampled_curr.csv");
    CSVHandler adaptivePaaHandler_volt("../log/PLAA_sampled_volt.csv");

// Random Sampling downsampling files
    CSVHandler randomHandler_soc("../log/random_sampled_soc.csv");
    CSVHandler randomHandler_curr("../log/random_sampled_curr.csv");
    CSVHandler randomHandler_volt("../log/random_sampled_volt.csv");

// RDP downsampling files
    CSVHandler RDPHandler_soc("../log/RDP_sampled_soc.csv");
    CSVHandler RDPHandler_curr("../log/RDP_sampled_curr.csv");
    CSVHandler RDPHandler_volt("../log/RDP_sampled_volt.csv");

    /*
// Wavelet downsampling files
    CSVHandler waveletHandler_soc("../log/DWT_sampled_soc.csv");
    CSVHandler waveletHandler_curr("../log/DWT_sampled_curr.csv");
    CSVHandler waveletHandler_volt("../log/DWT_sampled_volt.csv");
*/

// Save each downsampled dataset to its corresponding file with only two columns: time and the downsampled variable

// LTTB downsampled data
    lttbHandler_soc.saveToCSV(lttb_sampled_time_soc, lttb_sampled_soc);
    lttbHandler_curr.saveToCSV(lttb_sampled_time_curr, lttb_sampled_curr);
    lttbHandler_volt.saveToCSV(lttb_sampled_time_volt, lttb_sampled_volt);

// MinMax downsampled data
    minmaxHandler_soc.saveToCSV(MinMax_sampled_time_soc, MinMax_sampled_soc);
    minmaxHandler_curr.saveToCSV(MinMax_sampled_time_curr, MinMax_sampled_curr);
    minmaxHandler_volt.saveToCSV(MinMax_sampled_time_volt, MinMax_sampled_volt);

// PAA downsampled data
    paaHandler_soc.saveToCSV(PAA_sampled_time_soc, PAA_sampled_soc);
    paaHandler_curr.saveToCSV(PAA_sampled_time_curr, PAA_sampled_curr);
    paaHandler_volt.saveToCSV(PAA_sampled_time_volt, PAA_sampled_volt);

// Random Sampling downsampled data
    randomHandler_soc.saveToCSV(Random_sampled_time_soc, Random_sampled_soc);
    randomHandler_curr.saveToCSV(Random_sampled_time_curr, Random_sampled_curr);
    randomHandler_volt.saveToCSV(Random_sampled_time_volt, Random_sampled_volt);

// RDP downsampled data
    RDPHandler_soc.saveToCSV(RDP_sampled_time_soc, RDP_sampled_soc);
    RDPHandler_curr.saveToCSV(RDP_sampled_time_curr, RDP_sampled_curr);
    RDPHandler_volt.saveToCSV(RDP_sampled_time_volt, RDP_sampled_volt);

    /*
// Wavelet downsampled data
    waveletHandler_soc.saveToCSV(wavelet_sampled_time_soc, wavelet_sampled_soc);
    waveletHandler_curr.saveToCSV(wavelet_sampled_time_curr, wavelet_sampled_curr);
    waveletHandler_volt.saveToCSV(wavelet_sampled_time_volt, wavelet_sampled_volt);
*/
//#######################################################################################################

//################################ save the duration times ##############################################
    std::cout << std::endl;
    // Print the logfile with detailed durations for each variable
    std::ofstream logFile("/Users/carlocei/Desktop/DS_A_Tester/log/log.txt");

    logFile << "Largest Triangle Three Buckets duration (soc): " << lttb_time_soc << " seconds\n";
    logFile << "Largest Triangle Three Buckets duration (curr): " << lttb_time_curr << " seconds\n";
    logFile << "Largest Triangle Three Buckets duration (volt): " << lttb_time_volt << " seconds\n";

    logFile << "Min-Max Downsampling duration (soc): " << minmax_time_soc << " seconds\n";
    logFile << "Min-Max Downsampling duration (curr): " << minmax_time_curr << " seconds\n";
    logFile << "Min-Max Downsampling duration (volt): " << minmax_time_volt << " seconds\n";

    logFile << "Piecewise Aggregate Approximation duration (soc): " << paa_time_soc << " seconds\n";
    logFile << "Piecewise Aggregate Approximation duration (curr): " << paa_time_curr << " seconds\n";
    logFile << "Piecewise Aggregate Approximation duration (volt): " << paa_time_volt << " seconds\n";

    logFile << "Random Sampling duration (soc): " << random_sampling_time_soc << " seconds\n";
    logFile << "Random Sampling duration (curr): " << random_sampling_time_curr << " seconds\n";
    logFile << "Random Sampling duration (volt): " << random_sampling_time_volt << " seconds\n";

    logFile << "Douglas-Peucker (RDP) duration (soc): " << RDP_time_soc << " seconds\n";
    logFile << "Douglas-Peucker (RDP) duration (curr): " << RDP_time_curr << " seconds\n";
    logFile << "Douglas-Peucker (RDP) duration (volt): " << RDP_time_volt << " seconds\n";
/*
    logFile << "Discrete Wavelet Transform duration (soc): " << wavelet_time_soc << " seconds\n";
    logFile << "Discrete Wavelet Transform duration (curr): " << wavelet_time_curr << " seconds\n";
    logFile << "Discrete Wavelet Transform duration (volt): " << wavelet_time_volt << " seconds\n";
*/
    logFile.close();
//#######################################################################################################

//############################# launch the script to plot the graphs ####################################

    // Launch the Python script to plot the graphs
    std::string pythonScript = "/Users/carlocei/Desktop/DS_A_Tester/plotgraphs.py";  // Full path to script.py
    // Launch Python script with the path of the analyzed file as an argument
    std::string command = "python3 " + pythonScript + " " + selected_file;
    std::cout << command;
    int result = std::system(command.c_str());

    if (result != 0) {
        std::cerr << "Failed to execute Python script\n";
    }
//#######################################################################################################
    return 0;
}
