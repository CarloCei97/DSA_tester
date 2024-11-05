import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from datetime import datetime

# Ask the user whether they want to display the graphs or save them as screenshots
choice = input("\nDo you want to (1) display the graphs or (2) save screenshots? Enter 1 or 2: ")

# Check if a file path argument was provided
if len(sys.argv) > 1:
    original_file_path = sys.argv[1]
else:
    print("Error: No file path provided.")
    exit()

# Extract the original file name without the extension for naming the subfolder
original_file_name = os.path.splitext(os.path.basename(original_file_path))[0]

# Load the original CSV file
try:
    df1 = pd.read_csv(original_file_path)
    df1['time'] = pd.to_datetime(df1['time'])  # Convert 'time' to datetime
except pd.errors.EmptyDataError:
    print(f"Error: The file {original_file_path} is empty or invalid.")
    exit()

# Folder path for sampled data
sampled_folder = "/Users/carlocei/Desktop/DS_A_Tester/log/"  # Update with correct folder path for sampled signals

# Screenshot folder path based on the original file name
screenshot_folder = f"/Users/carlocei/Desktop/DS_A_Tester/screenshots/{original_file_name}"  # Update path as needed

# If saving screenshots, create the specific subfolder if it doesn't exist
if choice == '2':
    os.makedirs(screenshot_folder, exist_ok=True)

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
        formatted_title_original = original_file_name
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

        # Show or save based on user choice
        if choice == '1':
            plt.show()
        elif choice == '2':
            # Save the plot as a screenshot in the specific subfolder
            screenshot_filename = os.path.join(screenshot_folder, f"{sampled_file.replace('.csv', '')}_screenshot.png")
            plt.savefig(screenshot_filename, dpi=300)
            print(f"Saved screenshot as {screenshot_filename}")

        # Close the figure after saving or displaying
        plt.close(fig)

    except pd.errors.EmptyDataError:
        print(f"Error: The file {sampled_file} is empty or invalid.")
    except KeyError:
        print(f"Error: The file {sampled_file} is missing required columns. Skipping.")
