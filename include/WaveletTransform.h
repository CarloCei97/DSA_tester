#ifndef WAVELET_TRANSFORM_H
#define WAVELET_TRANSFORM_H

#include <iostream>
#include <vector>
#include <string>

class WaveletTransform {
public:
    // Method to downsample the time series and SoC while maintaining pairs
    static std::pair<std::vector<double>, std::vector<double>> downsample(
            const std::vector<double>& time_series,
            const std::vector<double>& SoC,
            const std::string& waveletName,
            double threshold
    );

private:
    // Haar Wavelet Transform (1 level of decomposition)
    static void haarWaveletTransform(const std::vector<double>& time_series, const std::vector<double>& SoC,
                              std::vector<double>& sampled_time_series, std::vector<double>& sampled_SoC);

    // Daubechies Wavelet Transform (1 level of decomposition)
    static void daubechiesWaveletTransform(const std::vector<double>& time_series, const std::vector<double>& SoC,
                                    std::vector<double>& sampled_time_series, std::vector<double>& sampled_SoC);

    // Apply thresholding to detail coefficients
    static void applyThreshold(std::vector<double>& detail, double threshold);
};

#endif // WAVELET_TRANSFORM_H
