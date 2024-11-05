import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from datetime import datetime

# Check if a file path argument was provided
if len(sys.argv) > 1:
    original_file_path = sys.argv[1]
else:
    print("Error: No file path provided.")
    exit()

# Load the original CSV file
try:
    df1 = pd.read_csv(original_file_path)
    df1['time'] = pd.to_datetime(df1['time'])  # Convert 'time' to datetime
except pd.errors.EmptyDataError:
    print(f"Error: The file {original_file_path} is empty or invalid.")
    exit()

# Folder path for sampled data
sampled_folder = "/Users/carlocei/Desktop/DS_A_Tester/log/"  # Update with correct folder path for sampled signals

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
        formatted_title_sampled = sampled_file.replace('_', ' ').replace('.csv', '')
        formatted_title_original = original_file_path.split('/')[-1].split('.')[0]
        formatted_title_dimension = formatted_title_sampled.split(' ')[-1]
        # Create subplots
        fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

        # First subplot for the original signal
        axs[0].plot(df1['time'], df1[formatted_title_dimension], marker=' ', linestyle='-', color='b')
        axs[0].set_title(f'(Original Signal - {formatted_title_original}  {formatted_title_dimension} )')
        axs[0].set_ylabel(formatted_title_dimension)
        axs[0].grid(True)

        # Second subplot for the sampled signal
        axs[1].plot(df2['Time'], df2['SoC'], marker=' ', linestyle='-', color='r')
        axs[1].set_title(f'(Sampled Signal - {formatted_title_sampled})')
        axs[1].set_xlabel('Time')
        axs[1].set_ylabel(formatted_title_dimension)
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
