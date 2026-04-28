import os
import requests
import re
import json

# Define the output folder
output_folder = r"C:\SynologyDrive\PhD\Projects\SpringBloom2021\Data\Gene_Data\Raw_Data\KEGG_Categories"
os.makedirs(output_folder, exist_ok=True)

# Define KEGG BRITE categories with subcategories
categories = {
    "Transporters": {"ko02000": "Transporters"},
    "Metabolism": {
        "ko01000": "Enzymes",
        "ko01001": "Protein kinases",
        "ko01009": "Protein phosphatases",
        "ko01002": "Peptidases and inhibitors",
        "ko01003": "Glycosyltransferases",
        "ko01005": "Lipopolysaccharide biosynthesis",
        "ko01011": "Peptidoglycan biosynthesis",
        "ko01004": "Lipid biosynthesis",
        "ko01008": "Polyketide biosynthesis",
        "ko01006": "Prenyltransferases",
        "ko01007": "Amino acid related enzymes",
        "ko00199": "Cytochrome P450",
        "ko00194": "Photosynthesis proteins"
    },
    "Genetic_Information_Processing": {
        "ko03000": "Transcription factors",
        "ko03021": "Transcription machinery",
        "ko03019": "Messenger RNA biogenesis",
        "ko03041": "Spliceosome",
        "ko03011": "Ribosome",
        "ko03009": "Ribosome biogenesis",
        "ko03016": "Transfer RNA biogenesis",
        "ko03012": "Translation factors",
        "ko03110": "Chaperones and folding catalysts",
        "ko04131": "Membrane trafficking",
        "ko04121": "Ubiquitin system",
        "ko03051": "Proteasome",
        "ko03032": "DNA replication proteins",
        "ko03036": "Chromosome proteins",
        "ko03400": "DNA repair proteins",
        "ko03029": "Mitochondrial biogenesis"
    },
    "Signaling_and_Cellular_Processes": {
        "ko02044": "Secretion system",
        "ko02042": "Bacterial toxins",
        "ko02022": "Two-component system",
        "ko02035": "Bacterial motility",
        "ko03037": "Cilium proteins",
        "ko04812": "Cytoskeleton proteins",
        "ko04147": "Exosome",
        "ko02048": "Prokaryotic defense",
        "ko04030": "G protein-coupled receptors",
        "ko04050": "Cytokine receptors",
        "ko04054": "Pattern recognition receptors",
        "ko03310": "Nuclear receptors",
        "ko04040": "Ion channels",
        "ko04031": "GTP-binding proteins",
        "ko04052": "Cytokines and neuropeptides",
        "ko04515": "Cell adhesion molecules",
        "ko04090": "CD molecules",
        "ko01504": "Antimicrobial resistance genes",
        "ko00535": "Proteoglycans",
        "ko00536": "Glycosaminoglycan binding proteins",
        "ko00537": "GPI-anchored proteins",
        "ko04091": "Lectins",
        "ko04990": "Domain-containing proteins"
    }
}

# Transporter subcategories
transporters_detailed = [
    "ABC transporters, prokaryotic type",
    "Major facilitator superfamily (MFS)",
    "Phosphotransferase system (PTS)",
    "Other transporters"
]

# Function to fetch KEGG BRITE data
def fetch_kegg_json(ko_id):
    url = f"https://rest.kegg.jp/get/br:{ko_id}/json"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"Failed to fetch {ko_id}: {response.status_code}")
            return None
    except Exception as e:
        print(f"Error fetching {ko_id}: {e}")
        return None

# Function to extract KO numbers from JSON
def extract_ko_numbers(data):
    ko_entries = []
    stack = [data]
    while stack:
        current = stack.pop()
        if "name" in current:
            match = re.search(r"\b(K\d{5})\b", current["name"])
            if match:
                ko_entries.append(match.group())
        if "children" in current:
            stack.extend(current["children"])
    return ko_entries

# Function to sanitize filenames
def sanitize_filename(name):
    return re.sub(r"[^\w\[\]-]", "_", name.replace(" ", "_"))

# Fetch and process each category
for category_name, subcategories in categories.items():
    print(f"Processing category: {category_name}")

    # Create category folder
    category_folder = os.path.join(output_folder, category_name)
    os.makedirs(category_folder, exist_ok=True)

    for ko_id, subcategory_name in subcategories.items():
        data = fetch_kegg_json(ko_id)
        if data:
            ko_numbers = extract_ko_numbers(data)
            output_file = os.path.join(category_folder, f"{subcategory_name}.txt")
            
            with open(output_file, "w") as f:
                for ko in ko_numbers:
                    f.write(f"{ko}\t{subcategory_name}\n")

            print(f"  - {subcategory_name}: {len(ko_numbers)} KO numbers saved.")

# Fetch and process transporters in more detail
print("\nProcessing detailed transporters...")

transporters_data = fetch_kegg_json("ko02000")  # Fetch full Transporters category
if transporters_data:
    for target_name in transporters_detailed:
        first_level_node = next(
            (child for child in transporters_data["children"] if child["name"] == target_name), None
        )
        if first_level_node:
            transporter_folder = os.path.join(output_folder, "Transporters", sanitize_filename(target_name))
            os.makedirs(transporter_folder, exist_ok=True)

            # Extract and save KO numbers
            for child in first_level_node.get("children", []):
                subcategory_name = child["name"]
                ko_numbers = extract_ko_numbers(child)
                output_file = os.path.join(transporter_folder, f"{sanitize_filename(subcategory_name)}.txt")
                
                with open(output_file, "w") as f:
                    for ko in ko_numbers:
                        f.write(f"{ko}\t{subcategory_name}\n")

                print(f"  - {subcategory_name}: {len(ko_numbers)} KO numbers saved.")
        else:
            print(f"⚠ No node found for '{target_name}'.")

print("\n✅ Processing complete!")
