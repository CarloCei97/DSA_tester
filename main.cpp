#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cstdlib> // Include for system() function
//#include <matplot/matplot.h>
#include "include/CSVHandler.h"
#include "include/DownsamplingAlgorithms.h"
#include "include/WaveletTransform.h"

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
    CSVHandler reader("/Users/carlocei/Desktop/DS_A_Tester/real_data/sample_data_SOC.csv");
    std::vector<CSVRecord> records = reader.readCSV();

    std::vector<double> time_series;
    std::vector<double> soc;

    for (const auto& record : records) {
        auto time_point = parseDateTime(record.datetime); // Convert string to time_point
        double time_value = timePointToDouble(time_point); // Convert time_point to double
        time_series.push_back(time_value);
        soc.push_back(record.soc_value);
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
    int result = system("python3 /Users/carlocei/Desktop/DS_A_Tester/plotgraphs.py");
    if (result != 0) {
        std::cerr << "Error: Failed to execute plotgraphs.py" << std::endl;
    }

    return 0;
}
