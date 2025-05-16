import pandas as pd
import matplotlib.pyplot as plt
import os

# Load data from a file
def load_data(file_path):
    """
    Load data from a tab-delimited file into a pandas DataFrame.

    Args:
        file_path (str): Path to the input file.

    Returns:
        pd.DataFrame: DataFrame containing the loaded data.
    """
    try:
        data = pd.read_csv(file_path, delimiter='\t')
        return data
    except Exception as e:
        print(f"Error loading data from {file_path}: {e}")
        return None

# Analyze and visualize data from multiple files
def analyze_data(files_and_labels):
    """
    Analyze and visualize the normalized data from multiple files.

    Args:
        files_and_labels (list of tuple): List of (DataFrame, label) pairs.
    """
    if not files_and_labels:
        print("No data to analyze.")
        return

    plt.figure(figsize=(10, 6))

    for data, label in files_and_labels:
        if data is None or data.empty:
            print(f"No data to plot for label: {label}")
            continue

        # Ensure required columns exist
        required_columns = ['w', f'{label}_exp', f'{label}_pred']
        if not all(col in data.columns for col in required_columns):
            print(f"Missing required columns in Dataset ({label}): {required_columns}")
            continue

        # Filter data within w range [80, 150]
        # data_filtered = data[(data['w'] >= 0) & (data['w'] <= 100)]
        data_filtered = data[data[f'{label}_exp'] > 0]

        if data_filtered.empty:
            print(f"No data points in range in the given range for label: {label}")
            continue

        # Normalize data[f'{label}_exp'] by its sum
        total_exp = data_filtered[f'{label}_exp'].sum()
        print(f'{label}:', total_exp)
        if total_exp == 0:
            print(f"Skipping {label} due to zero total sum.")
            continue

        normalized_values = data_filtered[f'{label}_exp'] / total_exp
        #normalized_values = data_filtered[f'{label}_exp']

        # Normalize data[f'{label}_pred'] by its sum
        total_pred = data_filtered[f'{label}_pred'].sum()
        if total_pred == 0:
            print(f"Skipping {label} due to zero total sum.")
            continue

        predicted_values = data_filtered[f'{label}_pred'] / total_pred
        # predicted_values = data_filtered[f'{label}_pred']
        
        # Scatter plot (dots) for experimental data
        plt.scatter(data_filtered['w'], normalized_values, label=f'{label}_exp (normalized)', s=30)

        # Line plot for predicted data (same color as corresponding scatter plot)
        plt.plot(data_filtered['w'], predicted_values, linestyle='-', alpha=0.8, label=f'{label}_pred')

    plt.xlabel('Hamming Weight (w)')
    plt.ylabel('Normalized Error Count')
    plt.yscale('log', base=2)
    plt.title('Experiment Results vs Predictions (w) (Normalized)')
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.show()

# Main function to load and analyze the data
if __name__ == "__main__":
    # Define the directory containing the data files
    output_dir = "../data/n1280_k640_c22"

    # List all relevant files in the directory
    files = [f for f in os.listdir(output_dir) if f.startswith("experiment_") and f.endswith(".txt")]

    if not files:
        print(f"No data files found in directory: {output_dir}")
    else:
        files_and_labels = []
        for file in files:
            file_path = os.path.join(output_dir, file)
            label = file.split("_")[1]  # Extract label from filename (e.g., "LB" from "experiment_LB_64.txt")
            data = load_data(file_path)
            if data is not None:
                files_and_labels.append((data, label))

        # Analyze and plot the data
        analyze_data(files_and_labels)
