import os
import csv
import re
import numpy as np
from post_processing_figure_4_2 import dimple_number  # Import the dimple_number function

# Initialize lists to store the extracted data
folder_names = []
width_list = []
height_list = []
alpha_list = []
ep_list = []
k_list = []
dimple_count_list = []  # List to store the number of dimples
dimple_area_list = []  # List to store the dimple area

# Get all folders in the current directory
folders = []
for width in np.linspace(10, 22, 13):
    for height in np.linspace(12, 43, 32):
        for per_alpha in np.linspace(0.15, 0.35, 9):
            alpha = np.arctan((height / 2) / (width * per_alpha))
            X = np.hstack([width, height, alpha])
            WorkerNumber = hash(tuple(X))

            folders.append(f"Files_{WorkerNumber}")

failed = 0

for folder in folders:
    folder_path = os.path.join(os.getcwd(), folder)

    if not os.path.exists(folder_path):
        print(f"Warning: Folder '{folder_path}' does not exist.")
        failed += 1
        continue

    ng_parameter_file = None
    depth_csv_file = None

    # Check for ngParameter python file and depth.csv file
    for file in os.listdir(folder_path):
        if "ngParameter" in file and file.endswith(".py"):
            ng_parameter_file = os.path.join(folder_path, file)
        if "depth.csv" in file:
            depth_csv_file = os.path.join(folder_path, file)

    # Default values
    width, height, alpha = None, None, None
    ep, k = 0, 0  # Default values for ep and k if depth.csv is not found
    dimple_count = 0  # Default value for dimple count
    dimple_area = 0  # Default value for dimple area

    # Extract variables from ngParameter python file
    if ng_parameter_file:
        with open(ng_parameter_file, 'r') as f:
            content = f.read()
            width_match = re.search(r'width\s*=\s*([\d.]+)', content)
            height_match = re.search(r'height\s*=\s*([\d.]+)', content)
            alpha_match = re.search(r'alpha\s*=\s*([\d.]+)', content)

            width = float(width_match.group(1)) if width_match else None
            height = float(height_match.group(1)) if height_match else None
            alpha = float(alpha_match.group(1)) if alpha_match else None
    else:
        print(f"Warning: No ngParameter file found in folder '{folder}'.")

    # Extract values from depth.csv file
    if depth_csv_file:
        with open(depth_csv_file, 'r') as f:
            csv_reader = csv.reader(f)
            next(csv_reader)  # Skip the header
            for row in csv_reader:
                ep = float(row[0])
                k = float(row[1])
    else:
        print(f"Warning: No depth.csv file found in folder '{folder}'. Setting ep and k to 0.")

    # Count the number of dimples and calculate dimple area
    if width and height and alpha:
        x_offset = height / (2 * np.tan(alpha))
        try:
            dimple_count, _, _ = dimple_number(folder_path, width, x_offset, verbose=False)
            dimple_area = width * height  # Calculate dimple area
        except Exception as e:
            print(f"Error: dimple_number function failed for folder '{folder}'. Error: {e}")
            failed += 1
            continue

    # Append the extracted values to the lists
    if ep != 0.1:
        failed += 1
    folder_names.append(folder)
    width_list.append(width)
    height_list.append(height)
    alpha_list.append(alpha)
    ep_list.append(ep)
    k_list.append(k)
    dimple_count_list.append(dimple_count)  # Append the dimple count
    dimple_area_list.append(dimple_area)  # Append the dimple area

print(f"Data extraction complete. {failed} folders failed.")

# Write the data to a CSV file
output_file = "extracted_data_4.csv"
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Folder", "Area", "Width", "Height", "Alpha", "Ep", "K", "N"])  # Header
    for i in range(len(folder_names)):
        writer.writerow([folder_names[i], dimple_area_list[i], width_list[i], height_list[i], np.rad2deg(alpha_list[i]), ep_list[i], k_list[i], dimple_count_list[i]])

print(f"Data extraction complete. Results saved to '{output_file}'.")