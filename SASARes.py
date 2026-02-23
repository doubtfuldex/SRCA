import freesasa
from Bio.PDB import PDBParser
import csv

# Define residue categories
hydrophobic = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
polar_uncharged = ['S', 'T', 'N', 'Q']
polar_charged = ['R', 'H', 'K', 'D', 'E']
glycine = ['G']  # Glycine

# Define a dictionary for converting three-letter to one-letter residue codes
three_to_one_map = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# Function to categorize residues
def categorize_residue(residue_code):
    if residue_code in hydrophobic:
        return 'Hydrophobic'
    elif residue_code in polar_uncharged:
        return 'Polar noncharged'
    elif residue_code in polar_charged:
        return 'Polar charged'
    elif residue_code in glycine:
        return 'Glycine'
    else:
        return 'Unknown'

# Load the structure from a PDB file
pdb_file = '1ubq.pdb'
structure = freesasa.Structure(pdb_file)

# Calculate the SASA
result = freesasa.calc(structure)

# Use Bio.PDB to parse residue names
parser = PDBParser(QUIET=True)
bio_structure = parser.get_structure('1UBQ', pdb_file)

# Get SASA for all residues
residue_sasa_dict = result.residueAreas()
sasa_dict_chain = residue_sasa_dict['A']  # Only chain A

# Prepare data for CSV
csv_data = [['Residue', 'SASA', 'Category']] 

# Collect residue data
for model in bio_structure:
    for chain in model:
        if chain.id == 'A': 
            for residue in chain:
                residue_id = residue.id[1] 
                residue_name = residue.get_resname() 
                one_letter_code = three_to_one_map.get(residue_name, '?') 

                # Retrieve SASA for the residue
                sasa_key = str(residue_id)  # Use residue number as string
                if sasa_key in sasa_dict_chain:
                    sasa = sasa_dict_chain[sasa_key].total
                    category = categorize_residue(one_letter_code)  # Determine category
                    # Add Residue name with number (e.g., M1), SASA, and category to CSV data
                    csv_data.append([f"{one_letter_code}{residue_id}", f"{sasa:.2f}", category])

# Write data to CSV file
output_file = 'residue_sasa.csv'
with open(output_file, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerows(csv_data)

print(f"Residue-wise SASA data has been written to {output_file}.")

