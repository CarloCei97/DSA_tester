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
std::pair<std::vector<double>, std::vector<double>> DownsamplingAlgorithms::largestTriangleThreeBuckets(const std::vector<double>& time_series,const std::vector<double>& soc,int n_out) {
    int n = soc.size(); // Get the size of the input data
    if (n_out >= n || n_out == 0) {
        // If the output size is greater than or equal to the input size or zero, return the original data
        return {time_series, soc};
    }

    std::vector<double> sampled_time(n_out); // Vector to store the downsampled time points
    std::vector<double> sampled_soc(n_out);  // Vector to store the downsampled SOC values

    sampled_time[0] = time_series[0]; // Set the first time point
    sampled_soc[0] = soc[0];          // Set the first SOC value
    sampled_time[n_out - 1] = time_series[n - 1]; // Set the last time point
    sampled_soc[n_out - 1] = soc[n - 1];          // Set the last SOC value

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
                double area = 0.5 * std::abs(
                        (time_series[j] - time_series[bucket_start]) * (soc[k] - soc[bucket_start]) -
                        (soc[j] - soc[bucket_start]) * (time_series[k] - time_series[bucket_start])
                );
                if (area > max_area) {
                    max_area = area;
                    best_idx = j;
                }
            }
        }

        sampled_time[i] = time_series[best_idx]; // Set the best time point for this bucket
        sampled_soc[i] = soc[best_idx];          // Set the corresponding SOC value for this bucket
    }

    return {sampled_time, sampled_soc}; // Return the downsampled time and SOC data
}

std::pair<std::vector<double>, std::vector<double>> DownsamplingAlgorithms::minmaxDownsampling(
        const std::vector<double>& time_series,
        const std::vector<double>& soc,
        int n_out
) {
    int n = soc.size(); // Get the size of the input data
    if (n_out >= n || n_out == 0) {
        // If the output size is greater than or equal to the input size or zero, return the original data
        return {time_series, soc};
    }

    std::vector<double> sampled_time; // Vector to store the downsampled time points
    std::vector<double> sampled_soc;  // Vector to store the downsampled SOC values

    // Calculate the bucket size based on the input and output size
    int bucket_size = n*5 / n_out;  //I added *5 cause the buckets were too small and the signal was trimmed.

    // Loop through each bucket to find the minimum and maximum values
    for (int i = 0; i < n; i += bucket_size) {
        auto bucket_start_time = time_series.begin() + i;
        auto bucket_end_time = (i + bucket_size < n) ? bucket_start_time + bucket_size : time_series.end();

        auto bucket_start_soc = soc.begin() + i;
        auto bucket_end_soc = (i + bucket_size < n) ? bucket_start_soc + bucket_size : soc.end();

        // Find the min and max SOC values and their corresponding time points
        auto min_it = std::min_element(bucket_start_soc, bucket_end_soc);
        auto max_it = std::max_element(bucket_start_soc, bucket_end_soc);

        // Calculate the indices for the min and max SOC values
        int min_index = std::distance(soc.begin(), min_it);
        int max_index = std::distance(soc.begin(), max_it);

        // Add the corresponding time and SOC values to the sampled vectors
        sampled_time.push_back(time_series[min_index]);
        sampled_soc.push_back(*min_it);

        sampled_time.push_back(time_series[max_index]);
        sampled_soc.push_back(*max_it);
    }

    // If the sampled data exceeds the output size, trim it to match n_out
    while (sampled_time.size() > n_out) {
        sampled_time.pop_back();
        sampled_soc.pop_back();
    }

    return {sampled_time, sampled_soc}; // Return the downsampled time and SOC data
}


// PAA algorithm
std::pair<std::vector<double>, std::vector<double>> DownsamplingAlgorithms::piecewiseAggregateApproximation(
        const std::vector<double>& time_series,
        const std::vector<double>& soc,
        int n_out
) {
    int n = soc.size(); // Get the size of the input data
    if (n_out >= n || n_out == 0) {
        // If the output size is greater than or equal to the input size or zero, return the original data
        return {time_series, soc};
    }

    std::vector<double> sampled_time(n_out); // Vector to store the downsampled time points
    std::vector<double> sampled_soc(n_out);  // Vector to store the downsampled SOC values

    int bucket_size = n / n_out; // Calculate the bucket size based on the input and output size

    // Loop through each bucket to calculate the average time point and SOC value
    for (int i = 0; i < n_out; ++i) {
        // Calculate average time in the bucket
        sampled_time[i] = std::accumulate(time_series.begin() + i * bucket_size,
                                          time_series.begin() + (i + 1) * bucket_size, 0.0) / bucket_size;

        // Calculate average SOC in the bucket
        sampled_soc[i] = std::accumulate(soc.begin() + i * bucket_size,
                                         soc.begin() + (i + 1) * bucket_size, 0.0) / bucket_size;
    }

    return {sampled_time, sampled_soc}; // Return the downsampled time and SOC data
}

//adaptive PAA
std::pair<std::vector<double>, std::vector<double>> DownsamplingAlgorithms::adaptivePAA(
        const std::vector<double>& time_series,
        const std::vector<double>& soc,
        int n_out,
        double variance_threshold
) {
    int n = soc.size(); // Get the size of the input data
    if (n_out >= n || n_out == 0) {
        // If the output size is greater than or equal to the input size or zero, return the original data
        return {time_series, soc};
    }

    std::vector<double> sampled_time; // Vector to store the downsampled time points
    std::vector<double> sampled_soc;  // Vector to store the downsampled SOC values
    int segment_length = n / n_out;   // Calculate the initial segment length based on the input and output size

    // Loop through the data and segment it based on variance
    for (int i = 0; i < n;) {
        // Get the current segment
        std::vector<double> segment_soc(soc.begin() + i, soc.begin() + std::min(i + segment_length, n));
        std::vector<double> segment_time(time_series.begin() + i, time_series.begin() + std::min(i + segment_length, n));

        // Calculate the mean and variance of the SOC values in the segment
        double variance = 0.0;
        double mean = std::accumulate(segment_soc.begin(), segment_soc.end(), 0.0) / segment_soc.size();
        for (double val : segment_soc) {
            variance += (val - mean) * (val - mean); // Calculate variance
        }
        variance /= segment_soc.size();

        // If the variance exceeds the threshold, break the segment into smaller pieces
        if (variance > variance_threshold) {
            int small_segment_length = segment_length / 2; // Use smaller segment length if variance is high
            std::vector<double> small_segment_soc(soc.begin() + i, soc.begin() + std::min(i + small_segment_length, n));
            std::vector<double> small_segment_time(time_series.begin() + i, time_series.begin() + std::min(i + small_segment_length, n));

            // Calculate the mean SOC value and the corresponding mean time value
            double small_segment_mean_soc = std::accumulate(small_segment_soc.begin(), small_segment_soc.end(), 0.0) / small_segment_soc.size();
            double small_segment_mean_time = std::accumulate(small_segment_time.begin(), small_segment_time.end(), 0.0) / small_segment_time.size();

            sampled_soc.push_back(small_segment_mean_soc);   // Add the mean SOC value of the small segment to the sampled SOC
            sampled_time.push_back(small_segment_mean_time); // Add the mean time of the small segment to the sampled time
            i += small_segment_length;
        } else {
            // Calculate the mean SOC and time for the larger segment
            double segment_mean_soc = std::accumulate(segment_soc.begin(), segment_soc.end(), 0.0) / segment_soc.size();
            double segment_mean_time = std::accumulate(segment_time.begin(), segment_time.end(), 0.0) / segment_time.size();

            sampled_soc.push_back(segment_mean_soc);   // Add the mean SOC value of the segment to the sampled SOC
            sampled_time.push_back(segment_mean_time); // Add the mean time of the segment to the sampled time
            i += segment_length;
        }

        if (sampled_soc.size() >= n_out) {
            break; // Stop if we have enough sampled data
        }
    }

    return {std::vector<double>(sampled_time.begin(), sampled_time.begin() + n_out),
            std::vector<double>(sampled_soc.begin(), sampled_soc.begin() + n_out)};
}

// random sampling algorithm
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

// hybrid PAA-MinMax
std::pair<std::vector<double>, std::vector<double>> DownsamplingAlgorithms::hybridPAAMinMax(
        const std::vector<double>& time_series,
        const std::vector<double>& soc,
        int n_out,
        double variance_threshold
) {
    int n = soc.size();
    if (n_out >= n || n_out == 0) {
        return {time_series, soc};
    }

    std::vector<double> sampled_time(n_out);
    std::vector<double> sampled_soc(n_out);
    int bucket_size = n / n_out;

    for (int i = 0; i < n_out; ++i) {
        // Extract the current segment for both SOC and time
        std::vector<double> segment_soc(soc.begin() + i * bucket_size, soc.begin() + std::min((i + 1) * bucket_size, n));
        std::vector<double> segment_time(time_series.begin() + i * bucket_size, time_series.begin() + std::min((i + 1) * bucket_size, n));

        // Calculate the average SOC value
        double avg_soc = std::accumulate(segment_soc.begin(), segment_soc.end(), 0.0) / segment_soc.size();

        // Check if the spike value exceeds the variance threshold
        double max_soc = *std::max_element(segment_soc.begin(), segment_soc.end());
        double spike_soc = (max_soc - avg_soc > variance_threshold) ? max_soc : avg_soc;

        // Find the corresponding time point for the max value if a spike is used
        double corresponding_time = (spike_soc == max_soc) ? segment_time[std::distance(segment_soc.begin(), std::max_element(segment_soc.begin(), segment_soc.end()))]
                                                           : std::accumulate(segment_time.begin(), segment_time.end(), 0.0) / segment_time.size();

        sampled_soc[i] = spike_soc;
        sampled_time[i] = corresponding_time;
    }

    return {sampled_time, sampled_soc};
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
