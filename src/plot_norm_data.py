import pandas as pd
import matplotlib.pyplot as plt
import os

# Load data from a file
def load_data(file_path):
    try:
        data = pd.read_csv(file_path, delimiter='\t')
        return data
    except Exception as e:
        print(f"Error loading data from {file_path}: {e}")
        return None

# Analyze, plot, and save normalized data
def analyze_and_save_data(files_and_labels, output_folder):
    if not files_and_labels:
        print("No data to analyze.")
        return

    os.makedirs(output_folder, exist_ok=True)
    plt.figure(figsize=(10, 6))

    for data, label in files_and_labels:
        if data is None or data.empty:
            print(f"No data for label: {label}")
            continue

        exp_col = f'{label}_exp'
        pred_col = f'{label}_pred'
        required_columns = ['w', exp_col, pred_col]
        if not all(col in data.columns for col in required_columns):
            print(f"Missing required columns in Dataset ({label}): {required_columns}")
            continue

        # data_filtered = data[data[exp_col] > 0]
        # data_filtered = data[data[pred_col] > 0.0001]
        data_filtered = data
        if data_filtered.empty:
            print(f"No valid data for label: {label}")
            continue

        total_exp = data_filtered[exp_col].sum()
        total_pred = data_filtered[pred_col].sum()

        if total_exp == 0 or total_pred == 0:
            print(f"Skipping {label} due to zero total.")
            continue

        # Overwrite original columns with normalized values
        data_filtered[exp_col] = data_filtered[exp_col] / total_exp
        data_filtered[pred_col] = data_filtered[pred_col] / total_pred

        # Plot
        plt.scatter(data_filtered['w'], data_filtered[exp_col],
                    label=f'{label}_exp (normalized)', s=30)
        plt.plot(data_filtered['w'], data_filtered[pred_col],
                 linestyle='-', alpha=0.8, label=f'{label}_pred')

        # Save normalized data using original column names
        output_data = data_filtered[['w', exp_col, pred_col]]
        output_file = os.path.join(output_folder, f'{label}.txt')
        output_data.to_csv(output_file, sep='\t', index=False)
        print(f"Saved normalized data for {label} to {output_file}")

    # Finalize plot
    plt.xlabel('Hamming Weight (w)')
    plt.ylabel('Normalized Error Count')
    plt.yscale('log', base=2)
    plt.title('Experiment Results vs Predictions (Normalized)')
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.show()

# Main execution
if __name__ == "__main__":
    input_dir = "../data/n1280_k640_c22"
    output_dir = "../data/n1280_k640_c22_normalized"

    files = [f for f in os.listdir(input_dir) if f.startswith("experiment_") and f.endswith(".txt")]

    if not files:
        print(f"No data files found in directory: {input_dir}")
    else:
        files_and_labels = []
        for file in files:
            file_path = os.path.join(input_dir, file)
            label = file.split("_")[1]  # e.g., "LB" from "experiment_LB_64.txt"
            data = load_data(file_path)
            if data is not None:
                files_and_labels.append((data, label))

        analyze_and_save_data(files_and_labels, output_dir)
