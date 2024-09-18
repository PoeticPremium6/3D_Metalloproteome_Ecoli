#############################
#Let's begin to annotate datasets.
# We will first add post-translational modifications to our binding data
#This was extracted from http://www.mutfunc.com/
import pandas as pd
#Load up datasets for annotation
binding_data = pd.read_csv(".../binding_data_New.csv")
non_metal_data = pd.read_csv(".../non_metal_data.csv")
ptm_data = pd.read_csv(".../MutFunc\\other_ptms.tab", sep='\t')  # Assuming it's tab-delimited

#Begin merging datasets with post-translational modifications
# Merge the PTM data with the binding data
common_ids_binding = set(binding_data['UniProt_ID']) & set(ptm_data['acc'])
common_ids_non_metal = set(non_metal_data['UniProt_ID']) & set(ptm_data['acc'])
# Extracting rows with common IDs
common_rows_binding = binding_data[binding_data['UniProt_ID'].isin(common_ids_binding)]
common_rows_non_metal = non_metal_data[non_metal_data['UniProt_ID'].isin(common_ids_non_metal)]
# Investigating the merged data with common IDs
merged_common_binding = pd.merge(common_rows_binding, ptm_data, left_on='UniProt_ID', right_on='acc', how='left')
merged_common_non_metal = pd.merge(common_rows_non_metal, ptm_data, left_on='UniProt_ID', right_on='acc', how='left')
#Make some descriptive statistics
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Count the modifications in both datasets
mod_count_binding = merged_common_binding['modification'].value_counts().reset_index()
mod_count_binding.columns = ['modification', 'Binding']

mod_count_non_binding = merged_common_non_metal['modification'].value_counts().reset_index()
mod_count_non_binding.columns = ['modification', 'Non-Binding']

# Merge the counts on modification
merged_counts = pd.merge(mod_count_binding, mod_count_non_binding, on='modification', how='outer').fillna(0)

# Sort the dataframe by Binding counts (you can change this to Non-Binding or use a combined score for custom sorting)
merged_counts = merged_counts.sort_values(by='Binding', ascending=False)

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
merged_counts.plot(x='modification', kind='bar', ax=ax, logy=True)  # logy=True sets the y-axis to a log scale
ax.set_title('Distribution of PTMs in Binding vs Non-Binding Data (Log Scale)')
ax.set_ylabel('Count (Log Scale)')
ax.set_xlabel('Modification Type')
# Rotate x-axis labels for better readability
plt.xticks(rotation=45, ha="right", rotation_mode="anchor")  # ha='right' aligns the end of the label instead of the center
# Ensure layout is tight
plt.tight_layout()
# Save the plot
save_path = '.../PTM_Distribution_Comparison.png'
try:
    plt.savefig(save_path, format='png', dpi=300)
    print(f"Plot saved successfully at {save_path}")
    plt.show()  # Only show the plot if save is successful
except Exception as e:
    print(f"Failed to save the plot due to: {str(e)}")


#Now let's load up our main datasets
#We should try to integrate MutFunc, Uniprot annotations with these data sets and prepare our mutation datasets
import pandas as pd
import csv

# Define the path to the transcription factor disruption data
tf_disruption_path = ".../MutFunc\\tfbs.tab"
# Define paths to the mutation data
ale_mutation_path = ".../ALE_Mutations.csv"
ltee_mutation_path = ".../LTEE_Mutations.csv"

# Function to parse mutation data
def parse_mutation_data(file_path):
    mutations = []
    with open(file_path, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter=';')
        for row in reader:
            mutations.append(row)
    return mutations

# Read the transcription factor disruption data
tf_disruption_df = pd.read_csv(tf_disruption_path, sep='\t')

# Parse ALE and LTEE mutations
ale_mutations = parse_mutation_data(ale_mutation_path)
ltee_mutations = parse_mutation_data(ltee_mutation_path)

# Convert the parsed mutation data to DataFrames
ale_mutation_df = pd.DataFrame(ale_mutations)
ltee_mutation_df = pd.DataFrame(ltee_mutations)

# Concatenate the mutation DataFrames
mutation_df = pd.concat([ale_mutation_df, ltee_mutation_df], ignore_index=True)

# Ensure the 'gene' and 'downstream' columns are in a comparable format
mutation_df['gene'] = mutation_df['gene'].str.strip().str.lower()
tf_disruption_df['downstream'] = tf_disruption_df['downstream'].str.strip().str.lower()

# Merge the DataFrames based on the 'gene' column in mutation_df and the 'downstream' column in tf_disruption_df
integrated_df = pd.merge(tf_disruption_df, mutation_df, left_on='downstream', right_on='gene')

# If there are no matching data points, print a message
if integrated_df.empty:
    print("No matching data found.")
else:
    print("Data has been successfully integrated!")
    print(integrated_df.head())

    # Specify the path to save the CSV file
    output_path = ".../integrated_data.csv"

    # Save the DataFrame as a CSV file
    integrated_df.to_csv(output_path, index=False)

    print(f"The merged data has been saved to {output_path}.")

#Let's collect Sequence data so we can annotate all proteins in Interpro

#Since Uniprot IDs are lacking from our mutation data, we need to fully annotated 'integrated_df' with these IDs
#So that it can map to both mutfunc data and our PDB database

import pandas as pd

# Load your original CSV file
input_file_path = ".../integrated_data.csv"
output_file_path = ".../blattner_numbers.csv"

# Read the CSV file
df = pd.read_csv(input_file_path)

# Extract the 'Blattner_Number' column
blattner_numbers = df['Blattner_Number']

# Save this column to a new CSV file
blattner_numbers.to_csv(output_file_path, index=False)

print(f"Blattner numbers have been saved to {output_file_path}.")

#It seems that the API in python can be quite slow, so we will do it on Uniprot website.
#Source database: Gene Name
#Target database: UniProtKB
#Now that we have the Uniprot ID numbers and desired info, we can now merge with varients
import pandas as pd

# Load the datasets
integrated_data_path = ".../integrated_data.csv"
idmapping_data_path = ".../idmapping.tsv"

integrated_df = pd.read_csv(integrated_data_path)
idmapping_df = pd.read_csv(idmapping_data_path, sep='\t')

# Rename the 'From' column in idmapping_df to 'Blattner_Number' for consistency
idmapping_df.rename(columns={'From': 'Blattner_Number'}, inplace=True)

# Merge the datasets on 'Blattner_Number'
master_df = pd.merge(integrated_df, idmapping_df, on='Blattner_Number', how='left')

# Save the merged dataset
output_path = ".../master_dataset.csv"
master_df.to_csv(output_path, index=False)

print("Merged dataset saved to", output_path)

#Let's start merging in mutfunc datasets again. We can start with 'other_ptms.tab'
import pandas as pd

# Load the datasets
master_dataset_path = ".../master_dataset.csv"
other_ptms_path = ".../MutFunc\\other_ptms.tab"

master_dataset = pd.read_csv(master_dataset_path)
other_ptms = pd.read_csv(other_ptms_path, sep="\t")  # other_ptms.tab is tab-separated

# Merge the datasets on UniProt ID
merged_dataset = pd.merge(master_dataset, other_ptms, left_on='Entry', right_on='acc', how='left')

# Save the merged dataset
output_path = ".../merged_dataset.csv"
merged_dataset.to_csv(output_path, index=False)

print("Merging completed. The merged dataset is saved at:", output_path)

#Alright! Let's clean this up
import pandas as pd

# Load the dataset
dataset_path = ".../merged_dataset.csv"
df = pd.read_csv(dataset_path)

# Remove unnecessary columns (edit this list based on your requirements)
columns_to_drop = ['entry', 'version', 'sequence_x',
                   'chr', 'pos', 'ref_strand', 'alt_strand',
                   'start', 'end', 'sequence_y',
                   'gene_strand', 'scored_strand',
                   'wt_score', 'mt_score','ic_difference','percentile_difference',
                   'cells', 'Reference Genome AA Seq']  # Replace with actual column names you want to drop
df.drop(columns_to_drop, axis=1, inplace=True)
print("Current Column Names:")
print(df.columns.tolist())

# Rename columns (if necessary)
df.rename(columns={'ref': 'TF_ref_Nucleotide', 'alt': 'TF_alt_Nucleotide',
                   'tf': 'Predicted_TF', 'downsteam': 'Gene_Downstream_TF',
                   'knockout_pvalue':'Knockout_pvalue', 'chip_evidence':'ChipSeq_PMID',
                   'score_difference':'Score_Difference_WildvsMutant', 'Blattner_Number': 'Locus_Tag',
                   'gene':'Gene_Name', 'AA_orig':'Evol_AA_orig', 'AA_residue':'Evol_AA_residue',
                   'AA_final': 'Evol_AA_final', 'Codon_Orig': 'Evol_Codon_Orig', 'Codon_Final': 'Evol_Codon_Final',
                   'del_ration':'Evol_Del_Ratio', 'Protein names': 'Protein_Name', 'Gene Names': 'Gene_Names',
                   'Length': 'Gene_Length', 'EC number': 'EC_Number', 'Interacts with': 'Known_PPI',
                   'Pathway':'Metabolic_Pathway', 'Gene Ontology IDs': 'GO_Terms',
                   'Subcellular location [CC]':'Subcellular_Location','PubMed ID':'PMID',
                   'Protein families':'Protein_Families', 'Domain [FT]': 'Domains',
                   'Topological domain':'Topo_Domains', 'acc':'Uniprot_ID',
                   'position':'Post-trans_Position', 'residue_type':'Post-trans_residue_type',
                   'pubmed':'Post-trans_PMID','modified_residue':'Post-trans_modified_residue',
                   'modification': 'Post-trans_modification'}, inplace=True)

# Reorder columns (if necessary)
df = df[['Uniprot_ID', 'Gene_Names', 'Gene_Name', 'Locus_Tag', 'Gene_Length', 'EC_Number', 'Entry', 'Protein_Name',
         'Known_PPI', 'Metabolic_Pathway', 'GO_Terms', 'Subcellular_Location',
         'PMID', 'Protein_Families', 'Domains', 'Topo_Domains',
         'TF_ref_Nucleotide', 'TF_alt_Nucleotide', 'impact', 'Predicted_TF', 'downstream',
         'Knockout_pvalue', 'ChipSeq_PMID', 'Score_Difference_WildvsMutant',
         'Evol_AA_orig', 'Evol_AA_residue', 'Evol_AA_final', 'Evol_Codon_Orig',
         'Evol_Codon_Final', 'del_ratio', 'Database', 'Post-trans_Position', 'Post-trans_residue_type',
         'Post-trans_PMID', 'Post-trans_modified_residue', 'Post-trans_modification']]

output_path = ".../cleaned_annotation_file.csv"
df.to_csv(output_path, index=False)
print("Cleaned file saved to", output_path)

