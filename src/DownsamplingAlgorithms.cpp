// DownsamplingAlgorithms.cpp

#include "../include/DownsamplingAlgorithms.h" // Include the corresponding header file
#include <cmath>    // Include cmath for mathematical functions like sin()
#include <algorithm> // Include algorithm for std::min_element, std::max_element, and std::shuffle
#include <cstdlib>  // Include cstdlib for functions like srand() and rand()

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

// Function to perform Largest Triangle Three Buckets (LTTB) downsampling
std::vector<double> DownsamplingAlgorithms::largestTriangleThreeBuckets(const std::vector<double>& data, int n_out) {
    int n = data.size(); // Get the size of the input data
    if (n_out >= n || n_out == 0) {
        // If the output size is greater than or equal to the input size or zero, return the original data
        return data;
    }

    std::vector<double> sampled(n_out); // Vector to store the downsampled data
    sampled[0] = data[0]; // Set the first element
    sampled[n_out - 1] = data[n - 1]; // Set the last element

    int bucket_size = (n - 2) / (n_out - 2); // Calculate the bucket size based on the input and output size

    // Loop through each bucket to find the best point
    for (int i = 1; i < n_out - 1; ++i) {
        int bucket_start = (i - 1) * bucket_size + 1;
        int bucket_end = i * bucket_size + 1;
        int next_bucket_start = bucket_end;
        int next_bucket_end = std::min((i + 1) * bucket_size + 1, n - 1);

        double max_area = -1;
        int best_idx = 0;
        // Find the point with the largest triangle area
        for (int j = bucket_start; j < bucket_end; ++j) {
            for (int k = next_bucket_start; k < next_bucket_end; ++k) {
                double area = 0.5 * std::abs(data[j] * (k - next_bucket_start) - data[k] * (j - bucket_start));
                if (area > max_area) {
                    max_area = area;
                    best_idx = j;
                }
            }
        }

        sampled[i] = data[best_idx]; // Set the best point for this bucket
    }

    return sampled; // Return the downsampled data
}

// Function to perform Min-Max downsampling
std::vector<double> DownsamplingAlgorithms::minmaxDownsampling(const std::vector<double>& data, int n_out) {
    int n = data.size(); // Get the size of the input data
    if (n_out >= n || n_out == 0) {
        // If the output size is greater than or equal to the input size or zero, return the original data
        return data;
    }

    std::vector<double> sampled; // Vector to store the downsampled data
    int bucket_size = n / n_out; // Calculate the bucket size based on the input and output size

    // Loop through each bucket to find the minimum and maximum values
    for (int i = 0; i < n; i += bucket_size) {
        auto bucket_start = data.begin() + i;
        auto bucket_end = (i + bucket_size < n) ? bucket_start + bucket_size : data.end();
        double min_value = *std::min_element(bucket_start, bucket_end);
        double max_value = *std::max_element(bucket_start, bucket_end);
        sampled.push_back(min_value); // Add minimum value to the sampled data
        sampled.push_back(max_value); // Add maximum value to the sampled data
    }

    // Return the downsampled data, limited to the output size
    return std::vector<double>(sampled.begin(), sampled.begin() + n_out);
}

// Function to perform Piecewise Aggregate Approximation (PAA) downsampling
std::vector<double> DownsamplingAlgorithms::piecewiseAggregateApproximation(const std::vector<double>& data, int n_out) {
    int n = data.size(); // Get the size of the input data
    if (n_out >= n || n_out == 0) {
        // If the output size is greater than or equal to the input size or zero, return the original data
        return data;
    }

    std::vector<double> sampled(n_out); // Vector to store the downsampled data
    int bucket_size = n / n_out; // Calculate the bucket size based on the input and output size

    // Loop through each bucket to calculate the average value
    for (int i = 0; i < n_out; ++i) {
        sampled[i] = std::accumulate(data.begin() + i * bucket_size, data.begin() + (i + 1) * bucket_size, 0.0) / bucket_size;
    }

    return sampled; // Return the downsampled data
}

// Function to perform Adaptive Piecewise Aggregate Approximation (PAA) downsampling
std::vector<double> DownsamplingAlgorithms::adaptivePAA(const std::vector<double>& data, int n_out, double variance_threshold) {
    int n = data.size(); // Get the size of the input data
    if (n_out >= n || n_out == 0) {
        // If the output size is greater than or equal to the input size or zero, return the original data
        return data;
    }

    std::vector<double> sampled; // Vector to store the downsampled data
    int segment_length = n / n_out; // Calculate the initial segment length based on the input and output size

    // Loop through the data and segment it based on variance
    for (int i = 0; i < n;) {
        std::vector<double> segment(data.begin() + i, data.begin() + std::min(i + segment_length, n));
        double variance = 0.0;
        double mean = std::accumulate(segment.begin(), segment.end(), 0.0) / segment.size();
        for (double val : segment) {
            variance += (val - mean) * (val - mean); // Calculate variance
        }
        variance /= segment.size();

        if (variance > variance_threshold) {
            int small_segment_length = segment_length / 2; // Use smaller segment length if variance is high
            double segment_mean = std::accumulate(data.begin() + i, data.begin() + std::min(i + small_segment_length, n), 0.0) / small_segment_length;
            sampled.push_back(segment_mean); // Add the mean of the small segment to the sampled data
            i += small_segment_length;
        } else {
            sampled.push_back(mean); // Add the mean of the larger segment to the sampled data
            i += segment_length;
        }

        if (sampled.size() >= n_out) {
            break; // Stop if we have enough sampled data
        }
    }

    return std::vector<double>(sampled.begin(), sampled.begin() + n_out); // Return the downsampled data
}

std::vector<double> DownsamplingAlgorithms::randomSampling(const std::vector<double>& data, int n_out) {
    int n = data.size();
    if (n_out >= n || n_out == 0) {
        return data;
    }

    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g);

    std::sort(indices.begin(), indices.begin() + n_out);
    std::vector<double> sampled(n_out);
    for (int i = 0; i < n_out; ++i) {
        sampled[i] = data[indices[i]];
    }

    return sampled;
}

std::vector<double> DownsamplingAlgorithms::hybridPAAMinMax(const std::vector<double>& data, int n_out, double variance_threshold) {
    int n = data.size();
    if (n_out >= n || n_out == 0) {
        return data;
    }

    std::vector<double> sampled(n_out);
    int bucket_size = n / n_out;

    for (int i = 0; i < n_out; ++i) {
        std::vector<double> segment(data.begin() + i * bucket_size, data.begin() + (i + 1) * bucket_size);
        double avg_value = std::accumulate(segment.begin(), segment.end(), 0.0) / segment.size();
        double spike_value = *std::max_element(segment.begin(), segment.end()) - avg_value > variance_threshold ? *std::max_element(segment.begin(), segment.end()) : avg_value;
        sampled[i] = spike_value;
    }

    return sampled;
}

// douglas algorithm
// Function to perform Douglas-Peucker downsampling
std::vector<double> DownsamplingAlgorithms::douglasPeucker(const std::vector<double>& data, int n_out, double epsilon) {
    int n = data.size(); // Get the size of the input data
    if (n_out >= n || n_out == 0) {
        // If the desired number of output points is greater than or equal to the input size
        // or if n_out is zero, return the original data
        return data;
    }

    std::vector<int> indexSet; // Vector to store the indices of the points that will be kept
    // Perform the recursive Douglas-Peucker algorithm to find which points to keep
    douglasPeuckerRecursive(data, 0, n - 1, epsilon, indexSet);

    // Create a vector to store the downsampled data
    std::vector<double> sampled;
    // Collect the data points corresponding to the indices in indexSet
    for (int idx : indexSet) {
        sampled.push_back(data[idx]);
    }

    return sampled; // Return the downsampled data
}

// Recursive function for Douglas-Peucker algorithm
void DownsamplingAlgorithms::douglasPeuckerRecursive(const std::vector<double>& data, int start, int end, double epsilon, std::vector<int>& indexSet) {
    if (end <= start) return; // Base case: If the segment is too small, do nothing

    // Initialize variables to find the point with the maximum perpendicular distance
    double maxDist = 0.0;
    int index = start;

    // Loop through the segment to find the point with the largest perpendicular distance from the line segment
    for (int i = start + 1; i < end; ++i) {
        double dist = perpendicularDistance(data, i, start, end); // Calculate the distance
        if (dist > maxDist) {
            maxDist = dist; // Update the maximum distance
            index = i; // Update the index of the point with the maximum distance
        }
    }

    // If the maximum distance is greater than the threshold epsilon
    if (maxDist > epsilon) {
        // Recursively apply the Douglas-Peucker algorithm to the segments on either side of the point with the maximum distance
        douglasPeuckerRecursive(data, start, index, epsilon, indexSet);
        indexSet.push_back(index); // Add the point with the maximum distance to the index set
        douglasPeuckerRecursive(data, index, end, epsilon, indexSet);
    }
}

// Utility function to calculate perpendicular distance
double DownsamplingAlgorithms::perpendicularDistance(const std::vector<double>& data, int index, int start, int end) {
    // Calculate the coefficients of the line equation (Ax + By + C = 0) for the line segment from start to end
    double A = data[end] - data[start];
    double B = start - end;
    double C = A * start + B * data[start];

    // Calculate and return the perpendicular distance from the point at 'index' to the line segment
    return std::abs(A * index + B * data[index] + C) / std::sqrt(A * A + B * B);
}

std::vector<double> DownsamplingAlgorithms::douglasPeuckerBinarySearch(const std::vector<double>& data, int n_out, double initialEpsilon) {
    double low = 0.0;
    double high = initialEpsilon;
    std::vector<double> sampled;
    int maxIterations = 100;  // Set a limit to prevent infinite looping
    int iteration = 0;
    double tolerance = 0.01 * n_out;  // Allow a small tolerance in the number of points

    while (iteration < maxIterations) {
        iteration++;
        double epsilon = (low + high) / 2.0;
        sampled = douglasPeucker(data, n_out, epsilon);

        if (sampled.size() > n_out + tolerance) {
            low = epsilon;  // Increase epsilon to reduce points
        } else if (sampled.size() < n_out - tolerance) {
            high = epsilon;  // Decrease epsilon to retain more points
        } else {
            break;  // If within tolerance, break the loop
        }
    }

    return sampled;
}
