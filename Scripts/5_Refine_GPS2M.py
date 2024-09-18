#Let's start to make some plots
#Let's start by making a smaller dataset to use for Figure 1 and Supplementary Figure 1
import pandas as pd
###METAL BINDING
import pandas as pd
def merge_datasets_in_chunks(metal_binding_path, master_annotation_path, chunksize=10000):
    columns_to_keep = [
        'UniProt_ID', 'Metal_type', 'Known_PPI', 'Metabolic_Pathway', 'GO_Terms',
        'Subcellular_Location', 'Protein_Families', 'Domains', 'Topo_Domains', 'EC_Number', 'Gene_Length'
    ]

    master_annotation_df = pd.read_csv(master_annotation_path, usecols=['Entry', 'Known_PPI', 'Metabolic_Pathway', 'GO_Terms',
                                                                        'Subcellular_Location', 'Protein_Families', 'Domains', 'Topo_Domains', 'EC_Number', 'Gene_Length'])
    # Keeping only the first occurrence of each UniProt ID
    master_annotation_df = master_annotation_df.drop_duplicates(subset='Entry')

    merged_chunks = []
    for chunk in pd.read_csv(metal_binding_path, chunksize=chunksize, usecols=['UniProt_ID', 'Metal_type']):
        merged_chunk = pd.merge(chunk, master_annotation_df, left_on='UniProt_ID', right_on='Entry', how='inner')
        merged_chunks.append(merged_chunk[columns_to_keep])

    return pd.concat(merged_chunks)

# Paths to your datasets
metal_binding_path = ".../binding_data_New.csv"
master_annotation_path = ".../categorized_dataset.csv"

# Merge datasets in chunks
refined_df = merge_datasets_in_chunks(metal_binding_path, master_annotation_path)

# Saving the output
output_path = ".../refined_integrated_dataset_metals_final.csv"
refined_df.to_csv(output_path, index=False)
print("Refined integrated dataset saved to", output_path)

####NON METAL BINDING
import pandas as pd

def merge_datasets_in_chunks(metal_binding_path, master_annotation_path, chunksize=10000):
    columns_to_keep = [
        'UniProt_ID', 'Known_PPI', 'Metabolic_Pathway', 'GO_Terms',
        'Subcellular_Location', 'Protein_Families', 'Domains', 'Topo_Domains', 'EC_Number', 'Gene_Length'
    ]

    master_annotation_df = pd.read_csv(master_annotation_path, usecols=['Entry', 'Known_PPI', 'Metabolic_Pathway', 'GO_Terms',
                                                                        'Subcellular_Location', 'Protein_Families', 'Domains', 'Topo_Domains', 'EC_Number', 'Gene_Length'])
    # Keeping only the first occurrence of each UniProt ID
    master_annotation_df = master_annotation_df.drop_duplicates(subset='Entry')

    merged_chunks = []
    for chunk in pd.read_csv(metal_binding_path, chunksize=chunksize, usecols=['UniProt_ID']):
        merged_chunk = pd.merge(chunk, master_annotation_df, left_on='UniProt_ID', right_on='Entry', how='inner')
        merged_chunks.append(merged_chunk[columns_to_keep])

    return pd.concat(merged_chunks)

# Paths to your datasets
metal_binding_path = ".../non_metal_data.csv"
master_annotation_path = ".../categorized_dataset.csv"

# Merge datasets in chunks
refined_df = merge_datasets_in_chunks(metal_binding_path, master_annotation_path)

# Saving the output
output_path = ".../refined_integrated_dataset_nonmetal.csv"
refined_df.to_csv(output_path, index=False)
print("Refined integrated dataset saved to", output_path)
