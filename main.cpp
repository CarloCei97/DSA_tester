#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cstdlib> // Include for system() function
//#include <matplot/matplot.h>
#include "include/CSVHandler.h"
#include "include/DownsamplingAlgorithms.h"

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
    const int n_out = 2000;
    const double variance_threshold = 5;
    const double epsilon = 0.5;

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

    // Perform downsampling using various methods and keeping track of the execution time
    auto start_time = std::chrono::high_resolution_clock::now();
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

    // Perform downsampling using Adaptive PAA and keeping track of the execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [APAA_sampled_time, APAA_sampled_soc] = DownsamplingAlgorithms::adaptivePAA(time_series, soc, n_out, variance_threshold);
    auto adaptive_paa_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    // Perform downsampling using Random Sampling and keeping track of the execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [Random_sampled_time, Random_sampled_soc] = DownsamplingAlgorithms::randomSampling(time_series, soc, n_out);
    auto random_sampling_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    // Perform downsampling using Hybrid PAA Min-Max and keeping track of the execution time
    start_time = std::chrono::high_resolution_clock::now();
    auto [hybrid_PAA_MinMax_sampled_time, hybrid_PAA_MinMax_sampled_soc] = DownsamplingAlgorithms::hybridPAAMinMax(time_series, soc, n_out, variance_threshold);
    auto hybrid_paa_minmax_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    auto Douglas_sampled = DownsamplingAlgorithms::douglasPeucker(soc, n_out, epsilon);
    auto Douglas_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();



    // Save the downsampled data to CSV files
    CSVHandler lttbHandler("../log/lttb_sampled.csv");
    CSVHandler minmaxHandler("../log/minmax_sampled.csv");
    CSVHandler paaHandler("../log/paa_sampled.csv");
    CSVHandler adaptivePaaHandler("../log/adaptive_paa_sampled.csv");
    CSVHandler randomHandler("../log/random_sampled.csv");
    CSVHandler hybridPaaMinMaxHandler("../log/hybrid_paa_minmax_sampled.csv");
    CSVHandler douglasHandler("../log/Douglas_sampled.csv");
    CSVHandler douglasAdaptiveHandler("../log/DouglasAdaptive_sampled.csv");
    lttbHandler.saveToCSV(lttb_sampled_time, lttb_sampled_soc);
    minmaxHandler.saveToCSV(MinMax_sampled_time, MinMax_sampled_soc);
    paaHandler.saveToCSV(PAA_sampled_time, PAA_sampled_soc);
    adaptivePaaHandler.saveToCSV(APAA_sampled_time, APAA_sampled_soc);
    randomHandler.saveToCSV(Random_sampled_time, Random_sampled_soc);
    hybridPaaMinMaxHandler.saveToCSV(hybrid_PAA_MinMax_sampled_time, hybrid_PAA_MinMax_sampled_soc);
    douglasHandler.saveToCSV(time_series, Douglas_sampled);

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
    std::ofstream logFile("/Users/carlocei/Desktop/NTU/Downsampling_algorithms_test/log/log.txt");
    logFile << "Largest Triangle Three Buckets duration: " << lttb_time << " seconds\n";
    logFile << "Min-Max Downsampling duration: " << minmax_time << " seconds\n";
    logFile << "Piecewise Aggregate Approximation duration: " << paa_time << " seconds\n";
    logFile << "Random Sampling duration: " << random_sampling_time << " seconds\n";
    logFile << "Adaptive PAA duration: " << adaptive_paa_time << " seconds\n";
    logFile << "Hybrid PAA MinMax duration: " << hybrid_paa_minmax_time << " seconds\n";
    logFile << "Douglas sampling duration: " << Douglas_time << " seconds\n";
    logFile.close();

    // Launch the Python script to plot the graphs
    int result = system("python3 /Users/carlocei/Desktop/DS_A_Tester/plotgraphs.py");
    if (result != 0) {
        std::cerr << "Error: Failed to execute plotgraphs.py" << std::endl;
    }

    return 0;
}
