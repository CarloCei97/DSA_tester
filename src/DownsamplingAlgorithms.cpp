// DownsamplingAlgorithms.cpp

#include "../include/DownsamplingAlgorithms.h" // Include the corresponding header file
#include <cmath>    // Include cmath for mathematical functions like sin()
#include <algorithm> // Include algorithm for std::min_element, std::max_element, and std::shuffle
#include <cstdlib>  // Include cstdlib for functions like srand() and rand()
#include <iostream>
#include <vector>
#include <cassert>
// Function to generate synthetic data for time series and voltage
void DownsamplingAlgorithms::generateSyntheticData(std::vector<double>& time_series, std::vector<double>& voltage) {
    unsigned seed = 42; // Set the seed for random number generation
    srand(seed); // Initialize the random number generator with the seed

    // Generate the time series data
    for (size_t i = 0; i < time_series.size(); ++i) {
        // Evenly distribute time points across the series from 0 to 1000
        time_series[i] = i * 1000.0 / (time_series.size() - 1);
    }

    // Generate synthetic voltage data with a combination of sine waves and random noise
    for (size_t i = 0; i < time_series.size(); ++i) {
        voltage[i] = sin(0.01 * time_series[i]) * 10 + sin(0.02 * time_series[i]) * 5 + ((rand() % 100) / 100.0 - 0.5);
    }

    // Introduce a downward spike in the voltage data at indices 2000 to 2020
    for (int i = 2000; i < 2020; ++i) {
        voltage[i] += (i - 2000) * (-10.0 / 20);
    }

    // Introduce an upward spike in the voltage data at indices 5000 to 5020
    for (int i = 5000; i < 5020; ++i) {
        voltage[i] += (i - 5000) * (10.0 / 20);
    }

    // Introduce a sharp drop in the voltage data at indices 7000 to 7010
    for (int i = 7000; i < 7010; ++i) {
        voltage[i] -= 5;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// LTTB ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to perform Largest Triangle Three Buckets (LTTB) downsampling
std::pair<std::vector<double>, std::vector<double>> DownsamplingAlgorithms::largestTriangleThreeBuckets(const std::vector<double>& time_series,const std::vector<double>& soc,int n_out) {
        size_t sourceSize = time_series.size();

        if (n_out >= sourceSize) {
            // If the requested output size is greater than or equal to the input size, return the original data
            return {time_series, soc};
        }

        if (n_out <= 2) {
            // If the requested output size is too small, just return the first and last points
            return { {time_series.front(), time_series.back()}, {soc.front(), soc.back()} };
        }

        std::vector<double> downsampled_time;
        std::vector<double> downsampled_soc;

        // Bucket size. Leave room for start and end data points
        double every = static_cast<double>(sourceSize - 2) / (n_out - 2);

        size_t aIndex = 0;  // Initially 'a' is the first point

        // Always add the first point
        downsampled_time.push_back(time_series[0]);
        downsampled_soc.push_back(soc[0]);

        for (size_t i = 0; i < n_out - 2; ++i) {
            // Calculate point average for the next bucket (containing 'c')
            double avgX = 0;
            double avgY = 0;
            size_t avgRangeStart = static_cast<size_t>((i + 1) * every) + 1;
            size_t avgRangeEnd = static_cast<size_t>((i + 2) * every) + 1;
            if (avgRangeEnd > sourceSize) {
                avgRangeEnd = sourceSize;
            }

            size_t avgRangeLength = avgRangeEnd - avgRangeStart;

            for (size_t j = avgRangeStart; j < avgRangeEnd; ++j) {
                avgX += time_series[j];
                avgY += soc[j];
            }
            avgX /= avgRangeLength;
            avgY /= avgRangeLength;

            // Get the range for this bucket
            size_t rangeOffs = static_cast<size_t>(i * every) + 1;
            size_t rangeTo = static_cast<size_t>((i + 1) * every) + 1;

            double pointAX = time_series[aIndex];
            double pointAY = soc[aIndex];

            double maxArea = -1;
            size_t nextAIndex = 0;

            // Find the point that forms the largest triangle with the current point 'a' and the averaged next bucket
            for (size_t j = rangeOffs; j < rangeTo; ++j) {
                double area = std::abs(
                        (pointAX - avgX) * (soc[j] - pointAY) -
                        (pointAX - time_series[j]) * (avgY - pointAY)
                );

                if (area > maxArea) {
                    maxArea = area;
                    nextAIndex = j;  // Next 'a' is this 'b'
                }
            }

            // Add the point forming the largest triangle
            downsampled_time.push_back(time_series[nextAIndex]);
            downsampled_soc.push_back(soc[nextAIndex]);

            // Update the index for the next 'a'
            aIndex = nextAIndex;
        }

        // Always add the last point
        downsampled_time.push_back(time_series.back());
        downsampled_soc.push_back(soc.back());

    return {downsampled_time, downsampled_soc};
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// MinMax ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::pair<std::vector<double>, std::vector<double>> DownsamplingAlgorithms::minmaxDownsampling(
        const std::vector<double>& time_series,
        const std::vector<double>& soc,
        int n_out
) {
    assert(n_out % 2 == 0);

    // If the output size is greater than or equal to the input size, return the original vectors
    if (n_out >= soc.size()) {
        return {time_series, soc};
    }

    std::vector<double> sampled_time;
    std::vector<double> sampled_soc;
    size_t block_size = (soc.size() - 1) / (n_out / 2);

    size_t start_idx = 0;
    for (int i = 0; i < n_out / 2; ++i) {
        size_t end_idx = std::min(start_idx + block_size, soc.size() - 1);

        auto [min_index, max_index] = argminmax(soc, start_idx, end_idx);

        if (min_index < max_index) {
            sampled_time.push_back(time_series[min_index]);
            sampled_soc.push_back(soc[min_index]);
            sampled_time.push_back(time_series[max_index]);
            sampled_soc.push_back(soc[max_index]);
        } else {
            sampled_time.push_back(time_series[max_index]);
            sampled_soc.push_back(soc[max_index]);
            sampled_time.push_back(time_series[min_index]);
            sampled_soc.push_back(soc[min_index]);
        }

        start_idx = end_idx;
    }

    return {sampled_time, sampled_soc};
}

std::pair<size_t, size_t> DownsamplingAlgorithms::argminmax(
        const std::vector<double>& arr,
        size_t start_idx,
        size_t end_idx
) {
    auto min_it = std::min_element(arr.begin() + start_idx, arr.begin() + end_idx);
    auto max_it = std::max_element(arr.begin() + start_idx, arr.begin() + end_idx);
    return {std::distance(arr.begin(), min_it), std::distance(arr.begin(), max_it)};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// PAA algorithm ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method for Piecewise Aggregate Approximation
std::pair<std::vector<double>, std::vector<double>> DownsamplingAlgorithms::piecewiseAggregateApproximation(
    const std::vector<double>& time_series,
    const std::vector<double>& soc,
    int n_out
) {
    int n = time_series.size();
    std::vector<double> reduced_time_series(n_out, 0.0);
    std::vector<double> reduced_soc(n_out, 0.0);
    int frameSize = n / n_out;

    // Calculate the mean of each frame for time series and soc
    // the mean of the time series equals the center
    for (int i = 0; i < n_out; ++i) {
        double sum_time_series = 0.0;
        double sum_soc = 0.0;

        for (int j = 0; j < frameSize; ++j) {
            sum_time_series += time_series[i * frameSize + j];
            sum_soc += soc[i * frameSize + j];
        }
        reduced_time_series[i] = sum_time_series / frameSize;
        reduced_soc[i] = sum_soc / frameSize;
    }

    return {reduced_time_series, reduced_soc};
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// PLAA ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::pair<std::vector<double>, std::vector<std::pair<double, double>>> DownsamplingAlgorithms::piecewiseLinearAggregateApproximation(
        const std::vector<double>& time_series,
        const std::vector<double>& soc,
        int n_out
) {
    int n = soc.size(); // Get the size of the input data
    if (n_out >= n || n_out == 0) {
        // If the output size is greater than or equal to the input size or zero, return the original data
        std::vector<std::pair<double, double>> soc_with_slopes;
        for (size_t i = 0; i < soc.size(); ++i) {
            soc_with_slopes.emplace_back(soc[i], 0.0); // Assume slope is 0 for original data
        }
        return {time_series, soc_with_slopes};
    }

    std::vector<double> sampled_time(n_out); // Vector to store the downsampled time points
    std::vector<std::pair<double, double>> sampled_soc(n_out);  // Vector to store the downsampled SOC values and slopes

    int bucket_size = n / n_out; // Calculate the bucket size based on the input and output size

    // Loop through each bucket to calculate the average time point, SOC value, and slope
    for (int i = 0; i < n_out; ++i) {
        // Calculate average time in the bucket
        sampled_time[i] = std::accumulate(time_series.begin() + i * bucket_size,
                                          time_series.begin() + (i + 1) * bucket_size, 0.0) / bucket_size;

        // Calculate average SOC in the bucket
        double mean_soc = std::accumulate(soc.begin() + i * bucket_size,
                                          soc.begin() + (i + 1) * bucket_size, 0.0) / bucket_size;

        // Calculate slope (linear regression) in the bucket
        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_xy = 0.0;
        double sum_x2 = 0.0;

        for (int j = 0; j < bucket_size; ++j) {
            int idx = i * bucket_size + j;
            double x = time_series[idx];
            double y = soc[idx];
            sum_x += x;
            sum_y += y;
            sum_xy += x * y;
            sum_x2 += x * x;
        }

        double slope = (bucket_size * sum_xy - sum_x * sum_y) / (bucket_size * sum_x2 - sum_x * sum_x);

        // Store the mean and slope for this bucket
        sampled_soc[i] = {mean_soc, slope};
    }

    return {sampled_time, sampled_soc}; // Return the downsampled time, SOC, and slopes
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// random sampling ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::pair<std::vector<double>, std::vector<double>> DownsamplingAlgorithms::randomSampling(
        const std::vector<double>& time_series,
        const std::vector<double>& soc,
        int n_out
) {
    int n = soc.size();
    if (n_out >= n || n_out == 0) {
        return {time_series, soc};
    }

    // Create a vector of indices
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0); // Fill the indices with values from 0 to n-1

    // Shuffle the indices randomly
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g);

    // Sort the first n_out indices to maintain the time order in the final result
    std::sort(indices.begin(), indices.begin() + n_out);

    // Create vectors to store the downsampled time and SOC values
    std::vector<double> sampled_time(n_out);
    std::vector<double> sampled_soc(n_out);

    // Populate the downsampled vectors using the selected random indices
    for (int i = 0; i < n_out; ++i) {
        sampled_time[i] = time_series[indices[i]];
        sampled_soc[i] = soc[indices[i]];
    }

    return {sampled_time, sampled_soc};
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Ramer-Douglas-Peucker ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::pair<std::vector<double>, std::vector<double>> DownsamplingAlgorithms::RDPsampling(
        const std::vector<double>& time_series,
        const std::vector<double>& soc,
        float epsilon
) {
    std::vector<double> time_result, soc_result;

    // Helper function to calculate perpendicular distance
    auto perpendicularDistance = [](double time_p, double soc_p, double time1, double soc1, double time2, double soc2) {
        double vec1_x = time_p - time1;
        double vec1_y = soc_p - soc1;
        double vec2_x = time2 - time1;
        double vec2_y = soc2 - soc1;
        double d_vec2 = sqrt(vec2_x * vec2_x + vec2_y * vec2_y);
        double cross_product = vec1_x * vec2_y - vec2_x * vec1_y;
        return fabs(cross_product / d_vec2);
    };

    // Find the point with the maximum distance
    double dmax = 0;
    int index = 0;
    for (int i = 1; i < time_series.size() - 1; ++i) {
        double d = perpendicularDistance(time_series[i], soc[i], time_series[0], soc[0], time_series[time_series.size() - 1], soc[soc.size() - 1]);
        if (d > dmax) {
            index = i;
            dmax = d;
        }
    }

    // If max distance is greater than epsilon, recursively simplify
    if (dmax > epsilon) {
        std::vector<double> pre_time_part, pre_soc_part, next_time_part, next_soc_part;

        for (int i = 0; i <= index; ++i) {
            pre_time_part.push_back(time_series[i]);
            pre_soc_part.push_back(soc[i]);
        }

        for (int i = index; i < time_series.size(); ++i) {
            next_time_part.push_back(time_series[i]);
            next_soc_part.push_back(soc[i]);
        }

        // Recursive calls
        auto result1 = RDPsampling(pre_time_part, pre_soc_part, epsilon);
        auto result2 = RDPsampling(next_time_part, next_soc_part, epsilon);

        // Combine results
        time_result.insert(time_result.end(), result1.first.begin(), result1.first.end());
        soc_result.insert(soc_result.end(), result1.second.begin(), result1.second.end());

        time_result.insert(time_result.end(), result2.first.begin() + 1, result2.first.end());
        soc_result.insert(soc_result.end(), result2.second.begin() + 1, result2.second.end());
    } else {
        // If not, just keep the first and last points
        time_result.push_back(time_series[0]);
        soc_result.push_back(soc[0]);
        time_result.push_back(time_series[time_series.size() - 1]);
        soc_result.push_back(soc[soc.size() - 1]);
    }

    return std::make_pair(time_result, soc_result);
}
