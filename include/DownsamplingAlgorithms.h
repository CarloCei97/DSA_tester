// DownsamplingAlgorithms.h
#ifndef DOWNSAMPLING_ALGORITHMS_H
#define DOWNSAMPLING_ALGORITHMS_H

#include <vector>
#include <numeric>
#include <random>

class DownsamplingAlgorithms {
public:
    static void generateSyntheticData(std::vector<double>& time_series, std::vector<double>& voltage);
    static std::vector<double> largestTriangleThreeBuckets(const std::vector<double>& data, int n_out);
    static std::vector<double> minmaxDownsampling(const std::vector<double>& data, int n_out);
    static std::vector<double> piecewiseAggregateApproximation(const std::vector<double>& data, int n_out);
    static std::vector<double> adaptivePAA(const std::vector<double>& data, int n_out, double variance_threshold);
    static std::vector<double> randomSampling(const std::vector<double>& data, int n_out);
    static std::vector<double> hybridPAAMinMax(const std::vector<double>& data, int n_out, double variance_threshold);
    // Douglas-Peucker algorithm method
    static std::vector<double> douglasPeucker(const std::vector<double>& data, int n_out, double epsilon);

private:
    // Helper functions for Douglas-Peucker
    static void douglasPeuckerRecursive(const std::vector<double>& data, int start, int end, double epsilon, std::vector<int>& indexSet);
    static double perpendicularDistance(const std::vector<double>& data, int index, int start, int end);
};

#endif // DOWNSAMPLING_ALGORITHMS_H
