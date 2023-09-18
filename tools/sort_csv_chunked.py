import csv
import os
import tempfile
import shutil

chunk_size = 400000  # Adjust as needed based on available memory

# Define the input and output file paths
input_file = r'./data/mobility/bs.csv'
output_file = r'./data/mobility/bs_sorted.csv'

# Define the primary and secondary columns to sort by (0-based indices)
primary_sort_column = 2  # Replace with the desired primary column index
secondary_sort_column = 3  # Replace with the desired secondary column index

# Create a temporary directory for splitting and sorting chunks
temp_dir = tempfile.mkdtemp()

try:
    # Split the large CSV file into smaller chunks
    chunk_files = []
    with open(input_file, 'r', newline='') as infile:
        reader = csv.reader(infile)
        header = next(reader)  # Save the header
        current_chunk = [header]
        for i, row in enumerate(reader):
            current_chunk.append(row)
            if i > 0 and i % chunk_size == 0:
                chunk_file = os.path.join(
                    temp_dir, f'chunk_{i // chunk_size}.csv')
                with open(chunk_file, 'w', newline='') as chunk:
                    writer = csv.writer(chunk)
                    writer.writerows(current_chunk)
                chunk_files.append(chunk_file)
                current_chunk = []
        if current_chunk:
            chunk_file = os.path.join(
                temp_dir, f'chunk_{i // chunk_size + 1}.csv')
            with open(chunk_file, 'w', newline='') as chunk:
                writer = csv.writer(chunk)
                writer.writerows(current_chunk)
            chunk_files.append(chunk_file)

    # Merge the sorted chunks into a single sorted CSV file
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)

        # Write the header only once
        writer.writerow(header)

        # Merge sort the chunks
        merged_data = []
        for chunk_file in chunk_files:
            with open(chunk_file, newline='') as chunk:
                data = list(csv.reader(chunk))
                merged_data.extend(data)

        # Define a custom sorting key function
        def custom_sort_key(row):
            return (row[primary_sort_column], row[secondary_sort_column])

        # Sort the merged data based on the custom sorting key
        merged_data.sort(key=custom_sort_key)

        # Write the sorted merged data to the output file
        writer.writerows(merged_data)

    print(
        f'Huge CSV file sorted by column {primary_sort_column} (primary) and {secondary_sort_column} (secondary) and saved as {output_file}')

finally:
    # Clean up temporary files and directories
    shutil.rmtree(temp_dir, ignore_errors=True)
