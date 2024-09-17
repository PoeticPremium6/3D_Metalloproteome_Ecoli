##################################CASE STUDY ###################
#Next, we want to establish a dataset to identify enzymes of industrial importance
#These proteins require substantial mutations and
# metal binding sites and mutations should be close to binding site
#Wecan assume that this is selecting for some sort of change that is affecting phenotypes
#and potentially metal-binding sites
#From this we can make two more figures:
#One on the prominence and diversity of metal-binding sites in carbon metabolism
#Second, on enzymes that would be prime targets for industry

import pandas as pd
#Let's filter out metal-binding sites
def filter_metal_binding_sites(input_csv_path, output_csv_path):
    # Load the dataset
    df = pd.read_csv(input_csv_path)

    # Filter out entries with a distance of more than 2 (keeping only distances of 2 or less)
    filtered_df = df[df["Distance"] <= 5]

    # Save the filtered dataset to a new CSV file
    filtered_df.to_csv(output_csv_path, index=False)
if __name__ == "__main__":
    input_csv_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv'  # Update this path to your actual CSV file path
    output_csv_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\filtered_binding_data.csv'  # Update this to your desired output file path

    filter_metal_binding_sites(input_csv_path, output_csv_path)
    print("Filtered dataset saved to:", output_csv_path)

#Now let's annotate metal-binding info with bIDs so they can be matched with ALE and LTEE
def merge_with_blattner_numbers_and_drop_unmatched(filtered_csv_path, id_mapping_path, output_csv_path):
    # Load the filtered dataset
    filtered_df = pd.read_csv(filtered_csv_path)

    # Load the ID mapping dataset
    id_mapping_df = pd.read_csv(id_mapping_path, sep='\t')

    # Merge the two datasets on the UniProt ID and Entry columns
    merged_df = pd.merge(filtered_df, id_mapping_df, left_on='UniProt_ID', right_on='Entry', how='left')

    # Select relevant columns and rename the Gene Names column to Blattner Numbers for clarity
    final_df = merged_df[
        ['PDB_name', 'UniProt_ID', 'Metal_type', 'Metal_coord', 'Residue_name', 'Residue_chain', 'Residue_number',
         'Residue_coord', 'Distance', 'Gene Names']]
    final_df.rename(columns={'Gene Names': 'Blattner Numbers'}, inplace=True)

    # Drop rows where Blattner Numbers are missing
    final_df = final_df.dropna(subset=['Blattner Numbers'])

    # Save the merged and filtered dataset to a new CSV file
    final_df.to_csv(output_csv_path, index=False)
if __name__ == "__main__":
    filtered_csv_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\filtered_binding_data.csv'  # Update this path
    id_mapping_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\idmapping.tsv'  # Update this path
    output_csv_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\merged_binding_data_with_blattner_filtered.csv'  # Update this to your desired output file path

    merge_with_blattner_numbers_and_drop_unmatched(filtered_csv_path, id_mapping_path, output_csv_path)
    print("Merged and filtered dataset saved to:", output_csv_path)

import pandas as pd
# Load the datasets
filtered_binding_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\merged_binding_data_with_blattner_filtered.csv'
ale_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\ALE_Mutations.csv'
ltee_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\LTEE_Mutations.csv'

binding_df = pd.read_csv(filtered_binding_data_path)
ale_df = pd.read_csv(ale_data_path, sep=';')
ltee_df = pd.read_csv(ltee_data_path, sep=';')

# Combine ALE and LTEE data for ease of processing
mutations_df = pd.concat([ale_df, ltee_df], ignore_index=True)

# Define a function to calculate proximity and count mutations
def count_mutations_close_to_binding_sites(binding_row, mutations_df, proximity_threshold=5):
    # Extract Blattner numbers from the binding data row
    blattner_numbers = binding_row['Blattner Numbers'].split()

    # Filter mutations relevant to the current protein based on Blattner numbers
    relevant_mutations = mutations_df[mutations_df['Blattner_Number'].isin(blattner_numbers)]

    # Count how many of these mutations occur within the proximity threshold of the binding site
    close_mutations_count = sum(
        abs(relevant_mutations['AA_residue'] - binding_row['Residue_number']) <= proximity_threshold
    )

    return close_mutations_count


# Apply the function to each row in the binding data to count mutations within proximity
binding_df['close_mutations_count'] = binding_df.apply(
    count_mutations_close_to_binding_sites, mutations_df=mutations_df, axis=1
)

# Filter out rows with 'close_mutations_count' equal to 0
filtered_binding_df = binding_df[binding_df['close_mutations_count'] > 0]

# Save the filtered dataset to a new CSV file
output_filtered_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\filtered_binding_data_with_mutation_counts.csv'
filtered_binding_df.to_csv(output_filtered_path, index=False)

print(f"Filtered dataset saved to {output_filtered_path}")

# Assuming 'gene' is a column in binding_df that specifies the gene name for each row
# Calculate the average close_mutations_count for each gene
average_mutation_counts = binding_df.groupby('Blattner Numbers')['close_mutations_count'].mean()

# Print out the results
print("Gene Name and Average Close Mutation Count:")
print(average_mutation_counts)

output_filtered_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Mutated_Metalloproteins.csv'
average_mutation_counts.to_csv(output_filtered_path, index=True)
#Now let's begin to merge our Mutated metalloproteins with annotations
def merge_datasets_in_chunks(metal_binding_path, master_annotation_path, chunksize=10000):
    # Update columns_to_keep to include Blattner Numbers and close_mutations_count
    columns_to_keep = [
        'UniProt_ID', 'Metal_type', 'Known_PPI', 'Metabolic_Pathway', 'GO_Terms',
        'Subcellular_Location', 'Protein_Families', 'Domains', 'Topo_Domains',
        'EC_Number', 'Gene_Length', 'Blattner Numbers', 'close_mutations_count'
    ]

    master_annotation_df = pd.read_csv(master_annotation_path, usecols=[
        'Entry', 'Known_PPI', 'Metabolic_Pathway', 'GO_Terms',
        'Subcellular_Location', 'Protein_Families', 'Domains', 'Topo_Domains',
        'EC_Number', 'Gene_Length'
    ])
    # Keeping only the first occurrence of each UniProt ID
    master_annotation_df = master_annotation_df.drop_duplicates(subset='Entry')

    merged_chunks = []
    # Include 'Blattner Numbers' and 'close_mutations_count' in usecols
    for chunk in pd.read_csv(metal_binding_path, chunksize=chunksize, usecols=[
        'UniProt_ID', 'Metal_type', 'Blattner Numbers', 'close_mutations_count'
    ]):
        merged_chunk = pd.merge(chunk, master_annotation_df, left_on='UniProt_ID', right_on='Entry', how='inner')
        # Ensure the final DataFrame includes the additional columns
        merged_chunks.append(merged_chunk[columns_to_keep])

    return pd.concat(merged_chunks)

# Continue with the merging process and EC Class integration as before
import pandas as pd

# Paths to your datasets
metal_binding_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\filtered_binding_data_with_mutation_counts.csv'
master_annotation_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\categorized_dataset.csv'
ec_class_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\EC_Class.csv'

# Merge the metal binding/mutation data with the master annotation dataset
refined_df = merge_datasets_in_chunks(metal_binding_path, master_annotation_path)

# Load the EC Class data
ec_class_df = pd.read_csv(ec_class_path)

# Merge the EC Class information based on EC_Number
# Ensure EC_Number columns are in a compatible format (string or float) before merging
refined_df = pd.merge(refined_df, ec_class_df, left_on='EC_Number', right_on='EC_Number', how='left')

# Define the path for the output CSV file
output_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\refined_binding_and_annotations.csv'

# Save the refined DataFrame to a CSV file
refined_df.to_csv(output_path, index=False)

print(f'Refined dataset saved to: {output_path}')

#We need to annotate the GO Terms for our analysis
import pandas as pd

# Load the main dataset
main_df_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\refined_binding_and_annotations.csv'
main_df = pd.read_csv(main_df_path)

# Load the GO terms annotations
go_terms_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\annotated_metals_go_terms.csv'
go_terms_df = pd.read_csv(go_terms_path)


# Function to map GO terms to their descriptions and categories
def map_go_terms(row, go_terms_df):
    descriptions = []
    categories = []
    go_terms = row['GO_Terms'].split('; ')  # Adjust split logic if necessary

    for go_term in go_terms:
        match = go_terms_df[go_terms_df['GO_Term'] == go_term.strip()]
        if not match.empty:
            descriptions.append(match['Description'].values[0])
            categories.append(match['Category'].values[0])

    # Combine descriptions and categories into strings
    row['GO_Descriptions'] = '; '.join(set(descriptions))
    row['GO_Categories'] = '; '.join(set(categories))
    return row


# Map GO terms for each row in the main DataFrame
main_df = main_df.apply(lambda row: map_go_terms(row, go_terms_df), axis=1)

# Define the path for the enriched output CSV file
enriched_output_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\enriched_binding_and_annotations_5.csv'

# Save the enriched DataFrame to a CSV file
main_df.to_csv(enriched_output_path, index=False)

print(f'Enriched dataset saved to: {enriched_output_path}')
