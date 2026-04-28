import os
import pandas as pd
import re

# Correct the directory path
directory_path = r'.'

# helper function to get the numbers out
def extract_number(SeqstatLine):
    if not SeqstatLine.strip():  # Check for empty or blank lines
        return "NA"
    match = re.search(r'\d[\d,]*\.?\d*', SeqstatLine)
    if match:
        return float(match.group().replace(',', ''))
    else:
        return "NA"

# Check if the directory exists
if not os.path.exists(directory_path):
    raise FileNotFoundError(f"Directory does not exist: {directory_path}")

# List all .txt files in the directory
txt_files = [f for f in os.listdir(directory_path) if f.endswith('.txt')]

# Mapping of column names to line indices
line_mapping = {
    "Number_Of_Contigs": 7,
    "Number_Of_Bases": 8,
    "Min_Contig_Length": 9,
    "Max_Contig_Length": 10,
    "Avg_Contig_Length": 11
}

# define columns
columns = ["Sample", "Data_Type", "Number_Of_Contigs", "Number_Of_Bases", "Min_Contig_Length", "Max_Contig_Length", "Avg_Contig_Length"]

# generate dataframe
df = pd.DataFrame(columns=columns)

# Iterate over .txt files and read each into a Pandas DataFrame
for filename in txt_files:
    file_path = os.path.join(directory_path, filename)
    row = []

    # Add sample name
    row.append(filename.split(".")[0])

    # Add raw or qf
    row.append(filename.split(".")[1])

    # Read file lines
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Dynamically add values based on the mapping
    for column, line_index in line_mapping.items():
        row.append(extract_number(lines[line_index]))

    # Append the row to the DataFrame
    df = pd.concat([df, pd.DataFrame([row], columns=columns)], ignore_index=True)

# File path to save the .txt file
output_file = os.path.join(directory_path, "Assembly_Report.txt")

# Write the DataFrame to a .txt file
df.to_csv(output_file, sep=' ', index=False)

print(f"File saved to: {output_file}")
