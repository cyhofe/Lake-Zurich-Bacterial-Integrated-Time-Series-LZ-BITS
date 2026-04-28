import os
import pandas as pd
import re

# Correct the directory path
directory_path = r'.'

# Helper function to extract numbers
def extract_number(line):
    if not line.strip():  # Check for empty or blank lines
        return "NA"
    match = re.search(r'\d[\d,]*\.?\d*', line)
    if match:
        return float(match.group().replace(',', ''))
    else:
        return "NA"

# Check if the directory exists
if not os.path.exists(directory_path):
    raise FileNotFoundError(f"Directory does not exist: {directory_path}")

# -------- Process Seqstat Files --------
# List all .seqstat.txt files in the directory
seqstat_files = [f for f in os.listdir(directory_path) if f.endswith('seqstat.txt')]

# Mapping of column names to line indices for seqstat
seqstat_line_mapping = {
    "Number_Of_Contigs": 7,
    "Number_Of_Bases": 8,
    "Min_Contig_Length": 9,
    "Max_Contig_Length": 10,
    "Avg_Contig_Length": 11
}

# Define columns for seqstat DataFrame
seqstat_columns = ["Sample", "Number_Of_Contigs", "Number_Of_Bases", "Min_Contig_Length", "Max_Contig_Length", "Avg_Contig_Length"]

# Generate seqstat DataFrame
seqstat_df = pd.DataFrame(columns=seqstat_columns)

# Process seqstat files
for filename in seqstat_files:
    file_path = os.path.join(directory_path, filename)
    row = []

    # Add sample name
    row.append(filename.split(".")[0])

    # Read file lines
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Dynamically add values based on the mapping
    for column, line_index in seqstat_line_mapping.items():
        row.append(extract_number(lines[line_index]))

    # Append the row to the seqstat DataFrame
    seqstat_df = pd.concat([seqstat_df, pd.DataFrame([row], columns=seqstat_columns)], ignore_index=True)

# -------- Process MappedReads Files --------
# List all .MappedReads.txt files in the directory
mapped_reads_files = [f for f in os.listdir(directory_path) if f.endswith('MappedReads.txt')]

# Define columns for MappedReads DataFrame
mapped_reads_columns = ["Sample", "Mapped_Reads", "Percent_Mapped_Reads"]

# Generate MappedReads DataFrame
mapped_reads_df = pd.DataFrame(columns=mapped_reads_columns)

# Process MappedReads files
for filename in mapped_reads_files:
    file_path = os.path.join(directory_path, filename)
    row = []

    # Add sample name
    row.append(filename.split(".")[0])

    # Read file lines
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract mapped reads
    mapped_reads_match = re.search(r'^(\d+)', lines[4])  # Fifth line contains mapped reads
    mapped_reads = int(mapped_reads_match.group(1)) if mapped_reads_match else "NA"
    row.append(mapped_reads)

    # Extract percent mapped
    percent_match = re.search(r'\((\d+)\.\d+%', lines[4])  # Fifth line contains percentage mapped
    percent_mapped = int(percent_match.group(1)) if percent_match else "NA"
    row.append(percent_mapped)

    # Append the row to the MappedReads DataFrame
    mapped_reads_df = pd.concat([mapped_reads_df, pd.DataFrame([row], columns=mapped_reads_columns)], ignore_index=True)

# -------- Process Total Changes Files --------
# List all .changes.total.txt files in the directory
changes_files = [f for f in os.listdir(directory_path) if f.endswith('.changes.total.txt')]

# Define columns for changes DataFrame
changes_columns = ["Sample", "Total_Changes"]

# Generate changes DataFrame
changes_df = pd.DataFrame(columns=changes_columns)

# Process changes files
for filename in changes_files:
    file_path = os.path.join(directory_path, filename)
    row = []

    # Add sample name
    row.append(filename.split(".")[0])

    # Read file lines
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract total changes
    changes_match = re.search(r'^(\d+)', lines[0])
    total_changes = int(changes_match.group(1)) if changes_match else "NA"
    row.append(total_changes)

    # Append the row to the changes DataFrame
    changes_df = pd.concat([changes_df, pd.DataFrame([row], columns=changes_columns)], ignore_index=True)

# -------- Combine DataFrames --------
# Merge all DataFrames on "Sample"
combined_df = seqstat_df.merge(mapped_reads_df, on="Sample", how="outer").merge(changes_df, on="Sample", how="outer")

# File path to save the combined DataFrame
output_file = os.path.join(directory_path, "Correction_3_Pilon_QualityReport.txt")

# Write the combined DataFrame to a .txt file
combined_df.to_csv(output_file, sep=' ', index=False)

print(f"File saved to: {output_file}")
