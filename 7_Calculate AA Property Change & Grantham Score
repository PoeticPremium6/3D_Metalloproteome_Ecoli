#AA Property Change and Grantham Score
import pandas as pd
import re

# File paths
annotation_file_path = "/Users/josspa/GPS-M/cleaned_annotation_scored.csv"
id_mapping_file_path = "/Users/josspa/GPS-M/Test/idmapping.tsv"
metal_binding_file_path = "/Users/josspa/GPS-M/binding_data_New.csv"

# Amino acid properties
aa_properties = {
    'A': {'General': 'Nonpolar', 'Properties': {'Non-polar', 'Aliphatic', 'Tiny', 'Small', 'Hydrophobic'}},
    'G': {'General': 'Nonpolar', 'Properties': {'Non-polar', 'Tiny', 'Small', 'Hydrophobic'}},
    'I': {'General': 'Nonpolar', 'Properties': {'Non-polar', 'Aliphatic', 'Small', 'Hydrophobic'}},
    'L': {'General': 'Nonpolar', 'Properties': {'Non-polar', 'Aliphatic', 'Small', 'Hydrophobic'}},
    'M': {'General': 'Nonpolar', 'Properties': {'Non-polar', 'Hydrophobic', 'Sulfur-containing'}},
    'F': {'General': 'Nonpolar', 'Properties': {'Non-polar', 'Aromatic', 'Hydrophobic'}},
    'P': {'General': 'Nonpolar', 'Properties': {'Non-polar', 'Small', 'Hydrophobic'}},
    'W': {'General': 'Nonpolar', 'Properties': {'Non-polar', 'Aromatic', 'Hydrophobic'}},
    'V': {'General': 'Nonpolar', 'Properties': {'Non-polar', 'Aliphatic', 'Small', 'Hydrophobic'}},
    'R': {'General': 'Positive', 'Properties': {'Polar', 'Charged', 'Basic', 'Hydrophilic'}},
    'K': {'General': 'Positive', 'Properties': {'Polar', 'Charged', 'Basic', 'Hydrophilic'}},
    'H': {'General': 'Positive', 'Properties': {'Polar', 'Aromatic', 'Charged', 'Basic', 'Hydrophilic'}},
    'D': {'General': 'Negative', 'Properties': {'Polar', 'Charged', 'Acidic', 'Hydrophilic', 'Small'}},
    'E': {'General': 'Negative', 'Properties': {'Polar', 'Charged', 'Acidic', 'Hydrophilic'}},
    'S': {'General': 'Polar', 'Properties': {'Polar', 'Tiny', 'Small', 'Hydrophilic'}},
    'T': {'General': 'Polar', 'Properties': {'Polar', 'Tiny', 'Small', 'Hydrophilic'}},
    'N': {'General': 'Polar', 'Properties': {'Polar', 'Hydrophilic', 'Amide', 'Small'}},
    'Q': {'General': 'Polar', 'Properties': {'Polar', 'Hydrophilic', 'Amide'}},
    'Y': {'General': 'Polar', 'Properties': {'Polar', 'Aromatic', 'Hydrophilic'}},
    'C': {'General': 'Polar', 'Properties': {'Polar', 'Tiny', 'Small', 'Hydrophobic', 'Sulfur-containing'}}
}

# Function to categorize based on Grantham score
def categorize_score(score):
    if pd.isna(score):
        return 'unknown'
    elif score <= 50:
        return 'conservative'
    elif score <= 100:
        return 'moderately conservative'
    elif score <= 150:
        return 'moderately radical'
    else:
        return 'radical'

# Load the annotated data
df_annotation = pd.read_csv(annotation_file_path, low_memory=False)

# Load the ID mapping data
df_id_mapping = pd.read_csv(id_mapping_file_path, delimiter='\t', low_memory=False)

# Function to extract Blattner numbers from 'Gene_Names' column
def extract_blattner_numbers(gene_names):
    return ' '.join(re.findall(r'\bb\d+\b', str(gene_names)))

# Function to calculate property change
def calculate_property_change(row):
    orig_property = aa_properties.get(row['Evol_AA_orig'], 'Unknown')
    final_property = aa_properties.get(row['Evol_AA_final'], 'Unknown')
    return f"{orig_property}_to_{final_property}"

# Apply the function to extract Blattner numbers
df_annotation['Blattner_Numbers'] = df_annotation['Gene_Names'].apply(extract_blattner_numbers)

# Apply the function to calculate property change
df_annotation['Property_Change'] = df_annotation.apply(calculate_property_change, axis=1)

# Create a dictionary to map Blattner numbers to Uniprot IDs
blattner_to_uniprot = pd.Series(df_id_mapping['Entry'].values, index=df_id_mapping['From']).to_dict()

# Map Blattner numbers to Uniprot IDs
df_annotation['Uniprot_ID'] = df_annotation['Blattner_Numbers'].apply(lambda x: blattner_to_uniprot.get(x.split()[0], 'Not Found') if x else 'Not Found')

# Filter rows that contain Blattner numbers and retain relevant columns
result_df = df_annotation[df_annotation['Blattner_Numbers'] != ''][['Blattner_Numbers', 'Evol_AA_orig', 'Evol_AA_residue', 'Evol_AA_final', 'Property_Change', 'Grantham_Score', 'Uniprot_ID']]

# Apply the function to categorize Grantham score
result_df['Grantham_Category'] = result_df['Grantham_Score'].apply(categorize_score)

# Load metal-binding data
df_metal_binding = pd.read_csv(metal_binding_file_path, low_memory=False)

# Rename 'UniProt_ID' to 'Uniprot_ID' in metal-binding data for consistency
df_metal_binding.rename(columns={'UniProt_ID': 'Uniprot_ID'}, inplace=True)

# Merge the dataframes on 'Uniprot_ID'
merged_df = pd.merge(result_df, df_metal_binding[['Uniprot_ID', 'Metal_type', 'Residue_number']], on='Uniprot_ID', how='left')

# Filter out rows where 'Metal_type' is blank
filtered_df = merged_df.dropna(subset=['Metal_type'])

# Display the first few rows of the filtered dataframe
print(filtered_df.head())

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming you have already loaded and merged your data into the 'filtered_df' dataframe

# Convert 'Evol_AA_residue' and 'Residue_number' to numeric values
filtered_df['Evol_AA_residue'] = pd.to_numeric(filtered_df['Evol_AA_residue'], errors='coerce')
filtered_df['Residue_number'] = pd.to_numeric(filtered_df['Residue_number'], errors='coerce')

# Calculate the distance between mutations and metal-binding sites
filtered_df['Distance'] = abs(filtered_df['Evol_AA_residue'] - filtered_df['Residue_number'])

# Drop rows with NaN in 'Distance', 'Evol_AA_residue', or 'Residue_number'
filtered_df.dropna(subset=['Distance', 'Evol_AA_residue', 'Residue_number'], inplace=True)

# Create and save filtered dataframes based on distance criteria
distance_criteria = [5, 10, 20]
for distance in distance_criteria:
    filtered_distance_df = filtered_df[filtered_df['Distance'] <= distance]
    filtered_distance_df.to_csv(f"/Users/josspa/GPS-M/processed_merged_data_distance_{distance}.csv", index=False)
