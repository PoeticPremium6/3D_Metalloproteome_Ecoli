#Now that we have all PDB/CIF files,
#we should process them to understand binding, geometry, and metal cofactor proximity
#while sorting out non-metalloproteins

import os
import pandas as pd
from biopandas.pdb import PandasPdb

def get_nearby_residues(metal_row, atom_df, cutoff_distance=10.0):
    distances = ((atom_df['x_coord'] - metal_row['x_coord']) ** 2 +
                 (atom_df['y_coord'] - metal_row['y_coord']) ** 2 +
                 (atom_df['z_coord'] - metal_row['z_coord']) ** 2) ** 0.5

    return atom_df[distances <= cutoff_distance]

def extract_uniprot_id_from_pdb(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('DBREF'):
                fields = line.split()
                if "UNP" in fields:
                    uniprot_index = fields.index("UNP")
                    if uniprot_index + 1 < len(fields):
                        return fields[uniprot_index + 1]
    return "Not_Found"

pdb_directory = ".../PDB_Total"
files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]

# List of metals by their residue names in PDB files. You may need to update or modify this list.
metals_residue_names = [
    'AG', 'AL', 'AU', 'BA', 'BE', 'BI', 'CA', 'CD', 'CE',
    'CO', 'CR', 'CS', 'CU', 'DY', 'ER', 'EU', 'FE', 'FR',
    'GA', 'GD', 'HF', 'HG', 'HO', 'IN', 'IR', 'K', 'LA', 'LI', 'LU',
    'MG', 'MN', 'MO', 'NA', 'NB', 'ND', 'NI', 'OS', 'PB',
    'PD', 'PR', 'PT', 'PU', 'RA', 'RB', 'RE', 'RH', 'RU', 'SB',
    'SC', 'SM', 'SN', 'SR', 'TA', 'TB', 'TI', 'TL', 'TM',
    'U', 'V', 'W', 'YB', 'ZN', 'ZR'
]

binding_data = []
non_metal_data = []

for file in files:
    file_path = os.path.join(pdb_directory, file)

    pdf = PandasPdb().read_pdb(file_path)
    adf = pdf.df['ATOM']
    hetdf = pdf.df['HETATM']

    uniprot_id = extract_uniprot_id_from_pdb(file_path)

    metal_df = hetdf[hetdf['residue_name'].isin(metals_residue_names)]

    if metal_df.empty:
        non_metal_data.append({
            'PDB_name': file.split('.')[0],
            'UniProt_ID': uniprot_id
        })
    else:
        for _, metal_row in metal_df.iterrows():
            nearby_residues = get_nearby_residues(metal_row, adf)
            for _, residue_row in nearby_residues.iterrows():
                binding_data.append({
                    'PDB_name': file.split('.')[0],
                    'UniProt_ID': uniprot_id,
                    'Metal_type': metal_row['residue_name'],
                    'Metal_coord': (metal_row['x_coord'], metal_row['y_coord'], metal_row['z_coord']),
                    'Residue_name': residue_row['residue_name'],
                    'Residue_chain': residue_row['chain_id'],
                    'Residue_number': residue_row['residue_number'],
                    'Residue_coord': (residue_row['x_coord'], residue_row['y_coord'], residue_row['z_coord']),
                    'Distance': ((metal_row['x_coord'] - residue_row['x_coord']) ** 2 +
                                 (metal_row['y_coord'] - residue_row['y_coord']) ** 2 +
                                 (metal_row['z_coord'] - residue_row['z_coord']) ** 2) ** 0.5
                })

binding_df = pd.DataFrame(binding_data)
binding_df.to_csv(".../binding_data_New.csv", index=False)

non_metal_df = pd.DataFrame(non_metal_data)
non_metal_df.to_csv(".../non_metal_data.csv", index=False)
