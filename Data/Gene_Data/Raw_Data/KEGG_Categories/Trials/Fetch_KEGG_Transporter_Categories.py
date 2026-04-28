import requests
import os
import re

# Fetch the JSON data from the KEGG REST API
url = "https://rest.kegg.jp/get/br:ko02000/json"
response = requests.get(url)
data = response.json()  # Parse the JSON data

# Function to extract KEGG IDs recursively
def extract_kegg_ids(node):
    kegg_ids = []
    if "name" in node and node["name"].startswith("K"):
        # Split the name by spaces and keep only the first part (the KEGG ID)
        kegg_ids.append(node["name"].split(" ")[0])
    if "children" in node:
        for child in node["children"]:
            kegg_ids.extend(extract_kegg_ids(child))
    return kegg_ids

# Function to sanitize filenames
def sanitize_filename(name):
    # Replace spaces with underscores
    sanitized_name = name.replace(" ", "_")
    # Remove invalid characters for filenames (e.g., :, /, \, etc.)
    sanitized_name = re.sub(r"[^\w\[\]-]", "", sanitized_name)
    return sanitized_name

# Function to process all subcategories and save KEGG IDs
def process_categories(root_node, output_folder):
    if "children" in root_node:
        for child in root_node["children"]:
            if "name" in child:
                subcategory_name = child["name"]

                # Extract KEGG IDs for this subcategory
                kegg_ids = extract_kegg_ids(child)

                # Log categories with no KEGG IDs
                if not kegg_ids:
                    print(f"Warning: No KEGG IDs found for category '{subcategory_name}'.")

                # Ensure the output folder exists
                os.makedirs(output_folder, exist_ok=True)

                # Sanitize filename
                sanitized_name = sanitize_filename(subcategory_name)
                filename = os.path.join(output_folder, f"{sanitized_name}.txt")

                # Save the KEGG IDs to the text file (only IDs)
                with open(filename, "w") as file:
                    if kegg_ids:
                        file.write("\n".join(kegg_ids))
                    else:
                        file.write("")  # Write nothing if no KEGG IDs

# Main categories to target
main_categories = [
    "ABC transporters, prokaryotic type",
    "Major facilitator superfamily (MFS)",
    "Phosphotransferase system (PTS)",
    "Other transporters",
]

# Specify output folder
output_folder = r"C:\SynologyDrive\PhD\Projects\SpringBloom2021\Data\Gene_Data\Raw_Data\KEGG_Categories\kegg_transporter_categories"

# Process each main category
for target_name in main_categories:
    first_level_node = None
    for child in data["children"]:
        if child["name"] == target_name:
            first_level_node = child
            break

    if first_level_node:
        # Generate files for each subcategory under the main category
        category_folder = os.path.join(output_folder, sanitize_filename(target_name))
        process_categories(first_level_node, category_folder)
        print(f"Files generated successfully for '{target_name}' in '{category_folder}' folder.")
    else:
        print(f"No node found with name '{target_name}'.")
