import os
import requests
import re
import json

# Define the base folder to save the files
output_folder = "C:\SynologyDrive\PhD\Projects\SpringBloom2021\Data\Gene_Data\Raw_Data\KEGG_Categories"
os.makedirs(output_folder, exist_ok=True)  # Create folder if it doesn't exist

# Define categories with their names and KO numbers
# categories adapted from: https://www.genome.jp/kegg/brite.html#gene, only took transporters as separate category, are allocated in signaling and cellular processes 
categories = {
    "Transporters": [
        "ko02000",  # Transporters
    ],
    "Metabolism": [
        "ko01000",  # Enzymes
        "ko01001",  # Protein kinases
        "ko01009",  # Protein phosphatases and associated proteins
        "ko01002",  # Peptidases and inhibitors
        "ko01003",  # Glycosyltransferases
        "ko01005",  # Lipopolysaccharide biosynthesis proteins
        "ko01011",  # Peptidoglycan biosynthesis and degradation proteins
        "ko01004",  # Lipid biosynthesis proteins
        "ko01008",  # Polyketide biosynthesis proteins
        "ko01006",  # Prenyltransferases
        "ko01007",  # Amino acid related enzymes
        "ko00199",  # Cytochrome P450
        "ko00194"   # Photosynthesis proteins
    ],
    "Genetic_Information_Processing": [
        "ko03000",  # Transcription factors
        "ko03021",  # Transcription machinery
        "ko03019",  # Messenger RNA biogenesis
        "ko03041",  # Spliceosome
        "ko03011",  # Ribosome
        "ko03009",  # Ribosome biogenesis
        "ko03016",  # Transfer RNA biogenesis
        "ko03012",  # Translation factors
        "ko03110",  # Chaperones and folding catalysts
        "ko04131",  # Membrane trafficking
        "ko04121",  # Ubiquitin system
        "ko03051",  # Proteasome
        "ko03032",  # DNA replication proteins
        "ko03036",  # Chromosome and associated proteins
        "ko03400",  # DNA repair and recombination proteins
        "ko03029"   # Mitochondrial biogenesis
    ],
    "Signaling_and_Cellular_Processes": [
        "ko02044",  # Secretion system
        "ko02042",  # Bacterial toxins
        "ko02022",  # Two-component system
        "ko02035",  # Bacterial motility proteins
        "ko03037",  # Cilium and associated proteins
        "ko04812",  # Cytoskeleton proteins
        "ko04147",  # Exosome
        "ko02048",  # Prokaryotic defense system
        "ko04030",  # G protein-coupled receptors
        "ko04050",  # Cytokine receptors
        "ko04054",  # Pattern recognition receptors
        "ko03310",  # Nuclear receptors
        "ko04040",  # Ion channels
        "ko04031",  # GTP-binding proteins
        "ko04052",  # Cytokines and neuropeptides
        "ko04515",  # Cell adhesion molecules
        "ko04090",  # CD molecules
        "ko01504",  # Antimicrobial resistance genes
        "ko00535",  # Proteoglycans
        "ko00536",  # Glycosaminoglycan binding proteins
        "ko00537",  # Glycosylphosphatidylinositol (GPI)-anchored proteins
        "ko04091",  # Lectins
        "ko04990"   # Domain-containing proteins not elsewhere classified
    ]
}

# Create a folder to store the output files
output_folder = "C:\SynologyDrive\PhD\Projects\SpringBloom2021\Data\Gene_Data\Raw_Data\KEGG_Categories\kegg_brite_categories"
os.makedirs(output_folder, exist_ok=True)

# Function to extract KO numbers from nested JSON
def extract_ko_numbers(data):
    ko_numbers = []
    stack = [data]  # Use a stack to avoid recursion
    while stack:
        current = stack.pop()  # Get the next item to process
        if "name" in current:
            match = re.search(r"\bK\d{5}\b", current["name"])
            if match:
                ko_numbers.append(match.group())
        if "children" in current:
            stack.extend(current["children"])  # Add children to the stack
    return ko_numbers

# Function to fetch JSON data and extract KO numbers
def fetch_ko_numbers_from_kegg(ko):
    url = f"https://rest.kegg.jp/get/br:{ko}/json"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            return extract_ko_numbers(data)
        else:
            print(f"Failed to fetch data for {ko}: {response.status_code}")
            return []
    except Exception as e:
        print(f"Error fetching data for {ko}: {e}")
        return []

# Process each category
for category_name, ko_list in categories.items():
    print(f"Processing category: {category_name}")
    category_ko_numbers = set()  # Use a set to avoid duplicates
    for ko in ko_list:
        ko_numbers = fetch_ko_numbers_from_kegg(ko)
        category_ko_numbers.update(ko_numbers)
    
    # Write KO numbers to a file
    output_file = os.path.join(output_folder, f"{category_name}.txt")
    with open(output_file, "w") as f:
        for ko_number in sorted(category_ko_numbers):
            f.write(f"{ko_number}\n")
    print(f"Saved {len(category_ko_numbers)} KO numbers to {output_file}")

print("\nProcessing complete!")