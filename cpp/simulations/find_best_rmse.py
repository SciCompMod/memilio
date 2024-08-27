import os
# Define the folder path
folder_path = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/cluster_results/rmse_grid/grid_search"

# Initialize a list to store the RMSE values and parameters
rmse_values = []
parameters = []

# Iterate over each file in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".txt"):
        file_path = os.path.join(folder_path, filename)
        with open(file_path, "r") as file:
            # Read all lines in the file
            lines = file.readlines()
            for line in lines:
                # Check if the line contains the RMSE value
                if "RMSE:" in line:
                    # Extract the RMSE value
                    rmse = float(line.split("RMSE:")[1].strip())
                    rmse_values.append(rmse)
                    # Extract the parameters
                    parameters.append(line.split("RMSE:")[0].strip())

# Sort the RMSE values in ascending order
sorted_indices = sorted(range(len(rmse_values)), key=lambda k: rmse_values[k])
rmse_values = [rmse_values[i] for i in sorted_indices]
parameters = [parameters[i] for i in sorted_indices]

# Print the 5 best RMSE values and parameters
print("Top 5 RMSE values and parameters:")
for i in range(20):
    print(f"Parameters: {parameters[i]}, RMSE: {rmse_values[i]}")