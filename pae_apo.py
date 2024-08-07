import json
import os
import numpy as np

# Directory where the JSON files are stored
directory = '004_json'

# List of JSON files in the specified directory
json_files = [
    'APO_004_fasta_scores_alphafold2_multimer_v3_model_5_seed_000.json',
    'APO_004_fasta_scores_alphafold2_multimer_v3_model_4_seed_000.json',
    'APO_004_fasta_scores_alphafold2_multimer_v3_model_3_seed_000.json',
    'APO_004_fasta_scores_alphafold2_multimer_v3_model_2_seed_000.json',
    'APO_004_fasta_scores_alphafold2_multimer_v3_model_1_seed_000.json',
]

# Dictionary to store the data from each JSON file
json_data = {}
pae_data = {}

# Print the current working directory
print(f"Current working directory: {os.getcwd()}")

# Load JSON data into a dictionary
for file_name in json_files:
    file_path = os.path.join(directory, file_name)
    if os.path.exists(file_path):
        try:
            with open(file_path, 'r') as file:
                data = json.load(file)
                json_data[file_name] = data
        except FileNotFoundError:
            print(f"File not found: {file_path}")
        except json.JSONDecodeError:
            print(f"Invalid JSON file: {file_path}")
    else:
        print(f"File not found: {file_path}")

# Pretty-print the JSON data for verification
for file_name, data in json_data.items():
    print(f"\nData from {file_name}:")
    print(json.dumps(data, indent=4))  # Pretty-printing the JSON data

# Using Numpy to calculate interface PAE
target_seq = ['MKAAVLTLAVLFLTGSQARHFWQQDEPPQSPWDRVKDLATVYVDVLKDSGRDYVSQFEGSALGKQLNLKLLDNWDSVTSTFSKLREQLGPVTQEFWDNLEKETEGLRQEMSKDLEEVKAKVQPYLDDFQKKWQEEMELYRQKVEPLRAELQEGARQKLHELQEKLSPLGEEMRDRARAHVDALRTHLAPYSDELRQRLAARLEALKENGGARLAEYHAKATEHLSTLSEKAKPALEDLRQGLLPVLESFKVSFLSALEEYTKKLNTQ']

def count_L(sequence):
    return len(sequence[0])

L = count_L(target_seq)
print(f"Size of APO-I-A protein: {L}")

def pae_calc(pae, L):
    """Calculate interface PAE given a PAE matrix and length of the target."""
    interface_pae = (np.mean(pae[:L, L:]) + np.mean(pae[L:, :L])) / 2
    return interface_pae

for file_name, data in json_data.items():
    if 'pae' in data and isinstance(data['pae'], list):
        pae_matrix = np.array(data['pae'])
        interface_pae = pae_calc(pae_matrix, L)
        print(f"Interface PAE for {file_name}: {interface_pae}")
    else:
        print(f"No 'pae' data found in {file_name}")
