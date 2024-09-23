
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# Load the first CSV file
file_path1 = "/Users/carlocei/Desktop/NTU/Downsampling_algorithms_test/real_data/sample_data_SOC.csv"  # Update this with the correct file path
df1 = pd.read_csv(file_path1)

# Convert the 'time' column in the first file to datetime format
df1['time'] = pd.to_datetime(df1['time'])

# Load the second CSV file
file_path2 = "/Users/carlocei/Desktop/NTU/Downsampling_algorithms_test/log/minmax_sampled.csv"  # Update this with the correct file path
df2 = pd.read_csv(file_path2)

# Convert the 'Time' column in the second file to datetime format
# Assuming 'Time' is a Unix timestamp, you can convert it like this
df2['Time'] = pd.to_datetime(df2['Time'], unit='s')

# Create subplots
fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# First subplot
axs[0].plot(df1['time'], df1['soc'], marker='o', linestyle='-', color='b')
axs[0].set_title('SOC vs. Time (File 1)')
axs[0].set_ylabel('SOC')
axs[0].grid(True)

# Second subplot
axs[1].plot(df2['Time'], df2['SoC'], marker='o', linestyle='-', color='r')
axs[1].set_title('SOC vs. Time (File 2)')
axs[1].set_xlabel('Time')
axs[1].set_ylabel('SOC')
axs[1].grid(True)

# Adjust layout
plt.xticks(rotation=45)
plt.tight_layout()

# Show the plot
plt.show()
