#include "../include/WaveletTransform.h"

// Haar Wavelet Transform (1 level of decomposition)
void WaveletTransform::haarWaveletTransform(const std::vector<double>& time_series, const std::vector<double>& SoC,
                                            std::vector<double>& sampled_time_series, std::vector<double>& sampled_SoC) {
    size_t len = SoC.size();
    for (size_t i = 0; i < len; i += 2) {
        double avg_SoC = (SoC[i] + SoC[i + 1]) / 2.0;
        double corresponding_time = (time_series[i] + time_series[i + 1]) / 2.0;

        // Store the downsampled pairs
        sampled_SoC.push_back(avg_SoC);
        sampled_time_series.push_back(corresponding_time);
    }
}

// Daubechies Wavelet Transform (1 level of decomposition)
void WaveletTransform::daubechiesWaveletTransform(const std::vector<double>& time_series, const std::vector<double>& SoC,
                                                  std::vector<double>& sampled_time_series, std::vector<double>& sampled_SoC) {
    // Using Daubechies-2 coefficients for simplicity
    static const double c0 = (1 + std::sqrt(3)) / (4 * std::sqrt(2));
    static const double c1 = (3 + std::sqrt(3)) / (4 * std::sqrt(2));
    static const double c2 = (3 - std::sqrt(3)) / (4 * std::sqrt(2));
    static const double c3 = (1 - std::sqrt(3)) / (4 * std::sqrt(2));

    size_t len = SoC.size();
    for (size_t i = 0; i + 3 < len; i += 2) {
        double avg_SoC = c0 * SoC[i] + c1 * SoC[i + 1] + c2 * SoC[i + 2] + c3 * SoC[i + 3];
        double corresponding_time = (time_series[i] + time_series[i + 1] + time_series[i + 2] + time_series[i + 3]) / 4.0;

        // Store the downsampled pairs
        sampled_SoC.push_back(avg_SoC);
        sampled_time_series.push_back(corresponding_time);
    }
}

// Apply thresholding to detail coefficients
void WaveletTransform::applyThreshold(std::vector<double>& detail, double threshold) {
    for (double& coeff : detail) {
        if (std::abs(coeff) < threshold) {
            coeff = 0.0;
        }
    }
}

// Method to downsample the time series and SoC
std::pair<std::vector<double>, std::vector<double>> WaveletTransform::downsample(
        const std::vector<double>& time_series,
        const std::vector<double>& SoC,
        const std::string& waveletName,
        double threshold
) {
    std::vector<double> sampled_time_series, sampled_SoC, detail;

    if (waveletName == "Haar") {
        haarWaveletTransform(time_series, SoC, sampled_time_series, sampled_SoC);
    }
    else if (waveletName == "Daubechies") {
        daubechiesWaveletTransform(time_series, SoC, sampled_time_series, sampled_SoC);
    } else {
        std::cerr << "Unsupported wavelet type!" << std::endl;
        return {time_series, SoC};  // Return original if wavelet is not supported
    }

    // Apply thresholding to detail coefficients (optional)
    applyThreshold(detail, threshold);

    return {sampled_time_series, sampled_SoC};  // Return downsampled time series and SoC pairs
}
