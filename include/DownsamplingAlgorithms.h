// DownsamplingAlgorithms.h
#ifndef DOWNSAMPLING_ALGORITHMS_H
#define DOWNSAMPLING_ALGORITHMS_H

#include <vector>
#include <numeric>
#include <random>

class DownsamplingAlgorithms {
public:
    static void generateSyntheticData(std::vector<double>& time_series, std::vector<double>& voltage);
    static std::pair<std::vector<double>, std::vector<double>> largestTriangleThreeBuckets(const std::vector<double>& time_series, const std::vector<double>& soc, int n_out);
    static std::pair<std::vector<double>, std::vector<double>> minmaxDownsampling(const std::vector<double>& time_series, const std::vector<double>& soc, int n_out);
    static std::pair<size_t, size_t> argminmax (const std::vector<double>& arr,size_t start_idx,size_t end_idx);
    static std::pair<std::vector<double>, std::vector<double>> piecewiseAggregateApproximation(const std::vector<double>& time_series,const std::vector<double>& soc,int n_out);
    static std::pair<std::vector<double>, std::vector<std::pair<double, double>>> piecewiseLinearAggregateApproximation(const std::vector<double>& time_series,const std::vector<double>& soc,int n_out);
    static std::pair<std::vector<double>, std::vector<double>> randomSampling(const std::vector<double>& time_series, const std::vector<double>& soc, int n_out);
    // Douglas-Peucker algorithm method
    static std::pair<std::vector<double>, std::vector<double>> RDPsampling(const std::vector<double>& time_series, const std::vector<double>& soc, float epsilon);

private:
};

#endif // DOWNSAMPLING_ALGORITHMS_H
