#6. Transcription Factor Disruption due to Mutations
import pandas as pd
import pandas as pd

# Define the path for the large dataset
merged_dataset_file = "/Users/josspa/GPS-M/Test/merged_dataset.csv"

# Load only the 'gene' and 'Entry' columns from the dataset
columns_to_load = ['gene', 'Entry']
merged_dataset = pd.read_csv(merged_dataset_file, usecols=columns_to_load, low_memory=False)

# Remove duplicate rows based on the 'gene' column
merged_dataset.drop_duplicates(subset='gene', inplace=True)

# Save this as a new file or continue processing
subset_file = "/Users/josspa/GPS-M/Test/merged_dataset_subset.csv"
merged_dataset.to_csv(subset_file, index=False)

# Optionally, display the first few rows of the dataset
print("Subset of Merged Dataset:\n", merged_dataset.head())
############
grantham_matrix = {
    ('A', 'A'): 0, ('A', 'R'): 112, ('A', 'N'): 111, ('A', 'D'): 126, ('A', 'C'): 195,
    ('A', 'Q'): 91, ('A', 'E'): 107, ('A', 'G'): 60, ('A', 'H'): 86, ('A', 'I'): 94,
    ('A', 'L'): 96, ('A', 'K'): 106, ('A', 'M'): 84, ('A', 'F'): 113, ('A', 'P'): 27,
    ('A', 'S'): 99, ('A', 'T'): 58, ('A', 'W'): 148, ('A', 'Y'): 112, ('A', 'V'): 64,
    ('R', 'A'): 112, ('R', 'R'): 0, ('R', 'N'): 86, ('R', 'D'): 96, ('R', 'C'): 180,
    ('R', 'Q'): 43, ('R', 'E'): 54, ('R', 'G'): 125, ('R', 'H'): 29, ('R', 'I'): 97,
    ('R', 'L'): 102, ('R', 'K'): 26, ('R', 'M'): 91, ('R', 'F'): 97, ('R', 'P'): 103,
    ('R', 'S'): 110, ('R', 'T'): 71, ('R', 'W'): 101, ('R', 'Y'): 77, ('R', 'V'): 96,
    ('N', 'A'): 111, ('N', 'R'): 86, ('N', 'N'): 0, ('N', 'D'): 23, ('N', 'C'): 139,
    ('N', 'Q'): 46, ('N', 'E'): 42, ('N', 'G'): 80, ('N', 'H'): 68, ('N', 'I'): 149,
    ('N', 'L'): 153, ('N', 'K'): 94, ('N', 'M'): 142, ('N', 'F'): 158, ('N', 'P'): 91,
    ('N', 'S'): 46, ('N', 'T'): 65, ('N', 'W'): 174, ('N', 'Y'): 143, ('N', 'V'): 133,
    ('D', 'A'): 126, ('D', 'R'): 96, ('D', 'N'): 23, ('D', 'D'): 0, ('D', 'C'): 154,
    ('D', 'Q'): 61, ('D', 'E'): 45, ('D', 'G'): 94, ('D', 'H'): 81, ('D', 'I'): 168,
    ('D', 'L'): 172, ('D', 'K'): 101, ('D', 'M'): 160, ('D', 'F'): 177, ('D', 'P'): 108,
    ('D', 'S'): 65, ('D', 'T'): 85, ('D', 'W'): 181, ('D', 'Y'): 160, ('D', 'V'): 152,
    ('C', 'A'): 195, ('C', 'R'): 180, ('C', 'N'): 139, ('C', 'D'): 154, ('C', 'C'): 0,
    ('C', 'Q'): 154, ('C', 'E'): 170, ('C', 'G'): 159, ('C', 'H'): 174, ('C', 'I'): 198,
    ('C', 'L'): 198, ('C', 'K'): 202, ('C', 'M'): 196, ('C', 'F'): 205, ('C', 'P'): 169,
    ('C', 'S'): 112, ('C', 'T'): 149, ('C', 'W'): 215, ('C', 'Y'): 194, ('C', 'V'): 192,
    ('Q', 'A'): 91, ('Q', 'R'): 43, ('Q', 'N'): 46, ('Q', 'D'): 61, ('Q', 'C'): 154,
    ('Q', 'Q'): 0, ('Q', 'E'): 29, ('Q', 'G'): 87, ('Q', 'H'): 24, ('Q', 'I'): 109,
    ('Q', 'L'): 113, ('Q', 'K'): 53, ('Q', 'M'): 101, ('Q', 'F'): 116, ('Q', 'P'): 76,
    ('Q', 'S'): 68, ('Q', 'T'): 42, ('Q', 'W'): 130, ('Q', 'Y'): 99, ('Q', 'V'): 96,
    ('E', 'A'): 107, ('E', 'R'): 54, ('E', 'N'): 42, ('E', 'D'): 45, ('E', 'C'): 170,
    ('E', 'Q'): 29, ('E', 'E'): 0, ('E', 'G'): 98, ('E', 'H'): 40, ('E', 'I'): 134,
    ('E', 'L'): 138, ('E', 'K'): 56, ('E', 'M'): 126, ('E', 'F'): 140, ('E', 'P'): 93,
    ('E', 'S'): 80, ('E', 'T'): 65, ('E', 'W'): 152, ('E', 'Y'): 122, ('E', 'V'): 121,
    ('G', 'A'): 60, ('G', 'R'): 125, ('G', 'N'): 80, ('G', 'D'): 94, ('G', 'C'): 159,
    ('G', 'Q'): 87, ('G', 'E'): 98, ('G', 'G'): 0, ('G', 'H'): 98, ('G', 'I'): 135,
    ('G', 'L'): 138, ('G', 'K'): 127, ('G', 'M'): 127, ('G', 'F'): 153, ('G', 'P'): 42,
    ('G', 'S'): 56, ('G', 'T'): 59, ('G', 'W'): 184, ('G', 'Y'): 147, ('G', 'V'): 109,
    ('H', 'A'): 86, ('H', 'R'): 29, ('H', 'N'): 68, ('H', 'D'): 81, ('H', 'C'): 174,
    ('H', 'Q'): 24, ('H', 'E'): 40, ('H', 'G'): 98, ('H', 'H'): 0, ('H', 'I'): 94,
    ('H', 'L'): 99, ('H', 'K'): 32, ('H', 'M'): 87, ('H', 'F'): 100, ('H', 'P'): 77,
    ('H', 'S'): 89, ('H', 'T'): 47, ('H', 'W'): 115, ('H', 'Y'): 83, ('H', 'V'): 84,
    ('I', 'A'): 94, ('I', 'R'): 97, ('I', 'N'): 149, ('I', 'D'): 168, ('I', 'C'): 198,
    ('I', 'Q'): 109, ('I', 'E'): 134, ('I', 'G'): 135, ('I', 'H'): 94, ('I', 'I'): 0,
    ('I', 'L'): 5, ('I', 'K'): 102, ('I', 'M'): 10, ('I', 'F'): 21, ('I', 'P'): 95,
    ('I', 'S'): 142, ('I', 'T'): 89, ('I', 'W'): 61, ('I', 'Y'): 33, ('I', 'V'): 29,
    ('L', 'A'): 96, ('L', 'R'): 102, ('L', 'N'): 153, ('L', 'D'): 172, ('L', 'C'): 198,
    ('L', 'Q'): 113, ('L', 'E'): 138, ('L', 'G'): 138, ('L', 'H'): 99, ('L', 'I'): 5,
    ('L', 'L'): 0, ('L', 'K'): 107, ('L', 'M'): 15, ('L', 'F'): 22, ('L', 'P'): 98,
    ('L', 'S'): 145, ('L', 'T'): 92, ('L', 'W'): 61, ('L', 'Y'): 36, ('L', 'V'): 32,
    ('K', 'A'): 106, ('K', 'R'): 26, ('K', 'N'): 94, ('K', 'D'): 101, ('K', 'C'): 202,
    ('K', 'Q'): 53, ('K', 'E'): 56, ('K', 'G'): 127, ('K', 'H'): 32, ('K', 'I'): 102,
    ('K', 'L'): 107, ('K', 'K'): 0, ('K', 'M'): 95, ('K', 'F'): 102, ('K', 'P'): 103,
    ('K', 'S'): 121, ('K', 'T'): 78, ('K', 'W'): 110, ('K', 'Y'): 85, ('K', 'V'): 97,
    ('M', 'A'): 84, ('M', 'R'): 91, ('M', 'N'): 142, ('M', 'D'): 160, ('M', 'C'): 196,
    ('M', 'Q'): 101, ('M', 'E'): 126, ('M', 'G'): 127, ('M', 'H'): 87, ('M', 'I'): 10,
    ('M', 'L'): 15, ('M', 'K'): 95, ('M', 'M'): 0, ('M', 'F'): 28, ('M', 'P'): 87,
    ('M', 'S'): 135, ('M', 'T'): 81, ('M', 'W'): 67, ('M', 'Y'): 36, ('M', 'V'): 21,
    ('F', 'A'): 113, ('F', 'R'): 97, ('F', 'N'): 158, ('F', 'D'): 177, ('F', 'C'): 205,
    ('F', 'Q'): 116, ('F', 'E'): 140, ('F', 'G'): 153, ('F', 'H'): 100, ('F', 'I'): 21,
    ('F', 'L'): 22, ('F', 'K'): 102, ('F', 'M'): 28, ('F', 'F'): 0, ('F', 'P'): 114,
    ('F', 'S'): 155, ('F', 'T'): 103, ('F', 'W'): 40, ('F', 'Y'): 22, ('F', 'V'): 50,
    ('P', 'A'): 27, ('P', 'R'): 103, ('P', 'N'): 91, ('P', 'D'): 108, ('P', 'C'): 169,
    ('P', 'Q'): 76, ('P', 'E'): 93, ('P', 'G'): 42, ('P', 'H'): 77, ('P', 'I'): 95,
    ('P', 'L'): 98, ('P', 'K'): 103, ('P', 'M'): 87, ('P', 'F'): 114, ('P', 'P'): 0,
    ('P', 'S'): 74, ('P', 'T'): 38, ('P', 'W'): 147, ('P', 'Y'): 110, ('P', 'V'): 68,
    ('S', 'A'): 99, ('S', 'R'): 110, ('S', 'N'): 46, ('S', 'D'): 65, ('S', 'C'): 112,
    ('S', 'Q'): 68, ('S', 'E'): 80, ('S', 'G'): 56, ('S', 'H'): 89, ('S', 'I'): 142,
    ('S', 'L'): 145, ('S', 'K'): 121, ('S', 'M'): 135, ('S', 'F'): 155, ('S', 'P'): 74,
    ('S', 'S'): 0, ('S', 'T'): 58, ('S', 'W'): 177, ('S', 'Y'): 144, ('S', 'V'): 124,
    ('T', 'A'): 58, ('T', 'R'): 71, ('T', 'N'): 65, ('T', 'D'): 85, ('T', 'C'): 149,
    ('T', 'Q'): 42, ('T', 'E'): 65, ('T', 'G'): 59, ('T', 'H'): 47, ('T', 'I'): 89,
    ('T', 'L'): 92, ('T', 'K'): 78, ('T', 'M'): 81, ('T', 'F'): 103, ('T', 'P'): 38,
    ('T', 'S'): 58, ('T', 'T'): 0, ('T', 'W'): 128, ('T', 'Y'): 92, ('T', 'V'): 69,
    ('W', 'A'): 148, ('W', 'R'): 101, ('W', 'N'): 174, ('W', 'D'): 181, ('W', 'C'): 215,
    ('W', 'Q'): 130, ('W', 'E'): 152, ('W', 'G'): 184, ('W', 'H'): 115, ('W', 'I'): 61,
    ('W', 'L'): 61, ('W', 'K'): 110, ('W', 'M'): 67, ('W', 'F'): 40, ('W', 'P'): 147,
    ('W', 'S'): 177, ('W', 'T'): 128, ('W', 'W'): 0, ('W', 'Y'): 37, ('W', 'V'): 88,
    ('Y', 'A'): 112, ('Y', 'R'): 77, ('Y', 'N'): 143, ('Y', 'D'): 160, ('Y', 'C'): 194,
    ('Y', 'Q'): 99, ('Y', 'E'): 122, ('Y', 'G'): 147, ('Y', 'H'): 83, ('Y', 'I'): 33,
    ('Y', 'L'): 36, ('Y', 'K'): 85, ('Y', 'M'): 36, ('Y', 'F'): 22, ('Y', 'P'): 110,
    ('Y', 'S'): 144, ('Y', 'T'): 92, ('Y', 'W'): 37, ('Y', 'Y'): 0, ('Y', 'V'): 55,
    ('V', 'A'): 64, ('V', 'R'): 96, ('V', 'N'): 133, ('V', 'D'): 152, ('V', 'C'): 192,
    ('V', 'Q'): 96, ('V', 'E'): 121, ('V', 'G'): 109, ('V', 'H'): 84, ('V', 'I'): 29,
    ('V', 'L'): 32, ('V', 'K'): 97, ('V', 'M'): 21, ('V', 'F'): 50, ('V', 'P'): 68,
    ('V', 'S'): 124, ('V', 'T'): 69, ('V', 'W'): 88, ('V', 'Y'): 55, ('V', 'V'): 0
}

import pandas as pd

# Function to categorize based on Grantham score
def categorize_score(score):
    if score is None:
        return 'unknown'
    elif score <= 50:
        return 'conservative'
    elif score <= 100:
        return 'moderately conservative'
    elif score <= 150:
        return 'moderately radical'
    else:
        return 'radical'

# Calculate Grantham score for a mutation
def get_grantham_score(aa_orig, aa_final):
    return grantham_matrix.get((aa_orig, aa_final), None)

# Define file paths
ale_mutations_file = "/Users/josspa/GPS-M/Variants/ALE_Mutations.csv"
ltee_mutations_file = "/Users/josspa/GPS-M/Variants/LTEE_Mutations.csv"
tfbs_file = "/Users/josspa/GPS-M/MutFunc/tfbs.tab"
mapping_file = "/Users/josspa/GPS-M/Test/merged_dataset_subset.csv"  # Path to the subset file

# Load the datasets
ale_mutations = pd.read_csv(ale_mutations_file, sep=';')
ltee_mutations = pd.read_csv(ltee_mutations_file, sep=';')
tfbs_data = pd.read_csv(tfbs_file, sep='\t')

# Merge the ALE and LTEE datasets
combined_mutations = pd.concat([ale_mutations, ltee_mutations])

# Remove rows where the 'gene' column is empty
combined_mutations = combined_mutations[combined_mutations['gene'].notna()]

# Calculate Grantham score and categorize mutations (ALE/LTEE data)
combined_mutations['Grantham_Score_Experimental'] = combined_mutations.apply(
    lambda row: get_grantham_score(row['AA_orig'], row['AA_final']), axis=1)
combined_mutations['Mutation_Category_Experimental'] = combined_mutations['Grantham_Score_Experimental'].apply(categorize_score)

# Calculate Grantham score and categorize mutations (TFBS data)
tfbs_data['Grantham_Score_Predicted'] = tfbs_data.apply(
    lambda row: get_grantham_score(row['ref'], row['alt']), axis=1)
tfbs_data['Mutation_Category_Predicted'] = tfbs_data['Grantham_Score_Predicted'].apply(categorize_score)

# Convert 'gene' to string in mutation data and 'downstream' in TFBS data
combined_mutations['gene'] = combined_mutations['gene'].astype(str)
tfbs_data['downstream'] = tfbs_data['downstream'].astype(str)

# Merge the combined mutation dataset with the TFBS dataset
merged_data = pd.merge(combined_mutations, tfbs_data, left_on='gene', right_on='downstream')

# Select specific columns from the merged dataset
selected_columns = [
    'gene', 'Database', 'AA_orig', 'AA_final', 'Grantham_Score_Experimental',
    'Mutation_Category_Experimental', 'ref', 'alt', 'impact', 'tf', 'downstream',
    'knockout_pvalue', 'Grantham_Score_Predicted', 'Mutation_Category_Predicted'
]

# Create a new DataFrame with only the selected columns
filtered_data = merged_data[selected_columns].copy()

# Load the UniProt mapping dataset and convert gene names to upper case
mapping_data = pd.read_csv(mapping_file)
mapping_data['gene'] = mapping_data['gene'].str.upper()

# Convert 'gene' names in filtered_data to upper case
filtered_data['gene'] = filtered_data['gene'].str.upper()

# Merge the filtered dataset with the UniProt mapping dataset
merged_with_uniprot = pd.merge(filtered_data, mapping_data, on='gene', how='left')

# Display the first few rows of the merged dataset with UniProt IDs
print("Merged Dataset with UniProt IDs:\n", merged_with_uniprot.head())

import pandas as pd

# Load the metal-binding dataset
metal_data_path = '/Users/josspa/GPS-M/refined_integrated_dataset_metals.csv'
metal_data = pd.read_csv(metal_data_path)

# Rename 'UniProt_ID' to 'Entry' in metal_data for consistent merging
metal_data.rename(columns={'UniProt_ID': 'Entry'}, inplace=True)

# Merge the metal data with the merged_with_uniprot dataset
enhanced_data = pd.merge(merged_with_uniprot, metal_data[['Entry', 'Metal_type']], on='Entry', how='left')

# Fill NaN valuxes in 'Metal_type' with 'N/A'
enhanced_data['Metal_type'].fillna('N/A', inplace=True)

# Print all column names in the enhanced_data DataFrame
print("Columns in the enhanced_data DataFrame:")
print(enhanced_data.columns.tolist())

# Define the columns to keep
columns_to_keep = ['gene', 'Mutation_Category_Experimental', 'impact', 'tf',
                   'downstream', 'Grantham_Score_Predicted', 'Mutation_Category_Predicted',
                   'Entry', 'Metal_type']

# Filter the DataFrame to only include the specified columns
filtered_enhanced_data = enhanced_data[columns_to_keep]

# Optionally, display the first few rows of the filtered DataFrame
print("Filtered Enhanced Dataset:\n", filtered_enhanced_data.head())

# Define the path for the output pickle file
output_pickle_path = '/Users/josspa/GPS-M/Test/filtered_enhanced_data.pkl'

# Save the filtered DataFrame to a pickle file
filtered_enhanced_data.to_pickle(output_pickle_path)
print(f"Filtered dataset saved as a pickle file at: {output_pickle_path}")
