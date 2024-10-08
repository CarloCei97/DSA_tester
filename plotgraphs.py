import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime

# Original file path
original_file_path = "/Users/carlocei/Desktop/DS_A_Tester/real_data/sample_data_SOC.csv"  # Update this with the correct file path

# Load the original CSV file
try:
    df1 = pd.read_csv(original_file_path)
    df1['time'] = pd.to_datetime(df1['time'])  # Convert 'time' to datetime
except pd.errors.EmptyDataError:
    print(f"Error: The file {original_file_path} is empty or invalid.")
    exit()

# Folder path for sampled data
sampled_folder = "/Users/carlocei/Desktop/DS_A_Tester/log/"  # Update this with the correct folder path for sampled signals

# List all .csv files in the sampled folder
sampled_files = [f for f in os.listdir(sampled_folder) if f.endswith('.csv')]

# Loop over the sampled files
for sampled_file in sampled_files:
    sampled_path = os.path.join(sampled_folder, sampled_file)
    
    try:
        # Load the sampled CSV file
        df2 = pd.read_csv(sampled_path)

        # Check if the file is empty
        if df2.empty:
            print(f"Warning: The file {sampled_file} is empty. Skipping.")
            continue

        # Convert the 'Time' column in the sampled file to datetime format (assuming Unix timestamp)
        df2['Time'] = pd.to_datetime(df2['Time'], unit='s')

        # Format the sampled file name: remove '.csv' and replace underscores with spaces
        formatted_title = sampled_file.replace('_', ' ').replace('.csv', '')

        # Create subplots
        fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

        # First subplot for the original signal
        axs[0].plot(df1['time'], df1['soc'], marker=' ', linestyle='-', color='b')
        axs[0].set_title(f'SOC vs. Time (Original Signal)')
        axs[0].set_ylabel('SOC')
        axs[0].grid(True)

        # Second subplot for the sampled signal
        axs[1].plot(df2['Time'], df2['SoC'], marker=' ', linestyle='-', color='r')
        axs[1].set_title(f'SOC vs. Time (Sampled Signal - {formatted_title})')
        axs[1].set_xlabel('Time')
        axs[1].set_ylabel('SOC')
        axs[1].grid(True)

        # Adjust layout
        plt.xticks(rotation=45)
        plt.tight_layout()

        # Show the plot
        plt.show()

    except pd.errors.EmptyDataError:
        print(f"Error: The file {sampled_file} is empty or invalid.")
    except KeyError:
        print(f"Error: The file {sampled_file} is missing required columns. Skipping.")
