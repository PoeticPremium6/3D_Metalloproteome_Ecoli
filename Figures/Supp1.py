#Let's make a plot to count how many times metals interact with residues
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Path to your CSV file
csv_file_path = ".../binding_data_New.csv"

# Read the CSV data into a DataFrame
binding_df = pd.read_csv(csv_file_path)

# List of amino acid residues
amino_acid_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# Extract the first 3 characters from the "Residue_name" column to get 3-letter residues
binding_df['Three_Letter_Residue'] = binding_df['Residue_name'].str[:3]

# Filter the DataFrame to include only amino acid residues
amino_acid_binding_df = binding_df[binding_df['Three_Letter_Residue'].isin(amino_acid_residues)]

# Count the occurrences of each amino acid residue
residue_counts = amino_acid_binding_df['Three_Letter_Residue'].value_counts()

# Create a bar plot using a single color (purple)
plt.figure(figsize=(12, 6))
sns.set_style("whitegrid")  # Set seaborn style
ax = sns.barplot(x=residue_counts.index, y=residue_counts.values, color="purple")  # Use a single color (purple)

#plt.title('Distribution of Amino Acid Residues Interacting with Metals', fontsize=16, fontweight='bold')
plt.xlabel('Amino Acid Residue', fontsize=14, fontweight='bold')
plt.ylabel('Count', fontsize=14, fontweight='bold')

# Setting tick labels to be bold
plt.xticks(rotation=45, fontsize=12, fontweight='bold')
plt.yticks(fontsize=10, fontweight='bold')

# Formatting y-axis tick labels to be more readable (remove scientific notation if present)
ax.get_yaxis().set_major_formatter(plt.ScalarFormatter())
ax.ticklabel_format(style='plain', axis='y')

# If you still see '1e6' on top or unclear scaling, this forces each tick to use integer formatting:
ax.set_yticklabels(['{:,}'.format(int(y)) for y in ax.get_yticks()])

plt.tight_layout()

# Save the figure with reduced whitespace
output_figure_path = ".../A_amino_acid_residue_distribution.png"
plt.savefig(output_figure_path, bbox_inches='tight', pad_inches=0.1)

print(f"Figure saved to: {output_figure_path}")

#Now let's make a plot to show the amount of times each metal interacts with a amino acid resiude
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy  # Ensure scipy is installed

# Path to your CSV file
csv_file_path = ".../binding_data_New.csv"

# Read the CSV data into a DataFrame
binding_df = pd.read_csv(csv_file_path)

# List of amino acid residues
amino_acid_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# Extract the first 3 characters from the "Residue_name" column to get 3-letter residues
binding_df['Three_Letter_Residue'] = binding_df['Residue_name'].str[:3]

# Filter the DataFrame to include only amino acid residues
amino_acid_binding_df = binding_df[binding_df['Three_Letter_Residue'].isin(amino_acid_residues)]

# Create a pivot table to count interactions between amino acids and metals
interaction_counts = amino_acid_binding_df.pivot_table(index='Three_Letter_Residue', columns='Metal_type', values='PDB_name', aggfunc='count', fill_value=0)

# Apply log transformation (adding a small constant to avoid log(0))
interaction_counts_log = np.log1p(interaction_counts)

# Use seaborn's clustermap function to create a heatmap with hierarchical clustering
g = sns.clustermap(interaction_counts_log, cmap='YlGnBu', figsize=(12, 10), annot=False,
                   cbar_kws={'label': 'Log Scale (1+x)'}, dendrogram_ratio=(.1, .1),
                   cbar_pos=(1.02, .2, .03, .4), linewidths=.5, linecolor='black')
g.ax_heatmap.set_xlabel('Metal Type', fontsize=14, fontweight='bold')
g.ax_heatmap.set_ylabel('3-Letter Amino Acid', fontsize=14, fontweight='bold')

# Enhance the tick labels for readability
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, fontsize=12, fontweight='bold')
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=12, fontweight='bold')

# Bold the color bar label
cbar = g.cax
cbar.set_title('Log Scale (1+x)', fontweight='bold', fontsize=10)
cbar.yaxis.label.set_size(10)
cbar.yaxis.label.set_weight('bold')

# Save the figure
output_figure_path = ".../B_amino_acid_metal_interaction_counts_clustermap.png"
plt.savefig(output_figure_path, bbox_inches='tight', pad_inches=0.1)


#Next, we should analyze protein domains
import pandas as pd
import re

# Function to extract and expand domain information
def expand_domain_info(df, metal_type_column='Metal_type'):
    expanded_rows = []
    for _, row in df.iterrows():
        metal_type = row[metal_type_column]
        domains = row['Domains']
        if pd.isna(domains) or domains.strip() == '':
            expanded_rows.append([metal_type, 'NA'])
        else:
            # Extract and split the domains based on '/note='
            domain_notes = re.findall(r'/note="([^"]+)"', domains)
            for note in domain_notes:
                expanded_rows.append([metal_type, note])
    return pd.DataFrame(expanded_rows, columns=[metal_type_column, 'Domains'])
# Function to process the metal-binding dataset
def process_metal_binding(file_path):
    df = pd.read_csv(file_path, usecols=['Metal_type', 'Domains'])
    return expand_domain_info(df)
# Function to process the non-metal binding dataset
def process_non_metal_binding(file_path):
    df = pd.read_csv(file_path, usecols=['Domains'])
    df['Metal_type'] = 'non-metal'  # Assign 'non-metal' to each row
    return expand_domain_info(df, metal_type_column='Metal_type')

# Process both datasets
metal_binding_data = process_metal_binding('.../refined_integrated_dataset_metals_final.csv')
non_metal_binding_data = process_non_metal_binding('.../refined_integrated_dataset_nonmetal.csv')

# Combine both datasets
combined_data = pd.concat([metal_binding_data, non_metal_binding_data])

# Save the combined dataset with a unique filename
combined_data.to_csv('.../metal_nonmetal_protein_domains.csv', index=False)

print("Data processing completed. The expanded data is saved to metal_nonmetal_protein_domains.csv")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the protein domain data
protein_domain_df = pd.read_csv('.../metal_nonmetal_protein_domains.csv')

# Replace missing values in 'Domains' with 'NA'
protein_domain_df['Domains'] = protein_domain_df['Domains'].fillna('NA')

# Aggregate the data by 'Metal_type' and 'Domains' and count the occurrences
aggregated_df = protein_domain_df.groupby(['Metal_type', 'Domains']).size().reset_index(name='Count')

# Filter out metals with no values in their columns
non_zero_metals = aggregated_df[aggregated_df['Count'] > 0]['Metal_type'].unique()
aggregated_df = aggregated_df[aggregated_df['Metal_type'].isin(non_zero_metals)]

# Create a pivot table for the heatmap
pivot_df = aggregated_df.pivot(index='Domains', columns='Metal_type', values='Count').fillna(0)

# Remove metal columns with all zeros
pivot_df = pivot_df.loc[:, (pivot_df != 0).any(axis=0)]

# Select the top domains based on their overall count to reduce the number of rows
# Sum counts across metals for each domain, sort, and select top N
top_domains_count = pivot_df.sum(axis=1).sort_values(ascending=False).head(20)  # Adjust N as needed
pivot_df = pivot_df.loc[top_domains_count.index]

# Apply log transformation to the counts
pivot_df_log = np.log10(pivot_df + 1)

# Plotting
plt.figure(figsize=(12, 8))
ax = sns.heatmap(pivot_df_log, annot=False, cmap="YlGnBu", linewidths=.5, linecolor='black', fmt="s")

# Adjusting visual elements for readability
plt.xticks(rotation=45, fontsize=12, fontweight='bold')  # X-axis labels
plt.yticks(rotation=0, fontsize=12, fontweight='bold')  # Y-axis labels
plt.xlabel('Metal Type', fontsize=14, fontweight='bold')  # X-axis title
plt.ylabel('Protein Domains', fontsize=14, fontweight='bold')  # Y-axis title
#plt.title('Log-Scale Heatmap of Top Protein Domain Abundance by Metal Type', fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig('.../C_Top_Protein_Domains.png', dpi=300, bbox_inches='tight')
plt.show()

#Let's analyze metal-binding proteins and their subcellular location
import pandas as pd
import os
# Function to clean the 'Subcellular_Location' column
def clean_location(location):
    if pd.isna(location):
        return 'Unknown'
    if 'SUBCELLULAR LOCATION:' in location:
        location = location.split('SUBCELLULAR LOCATION:')[1].split('.')[0].strip()
        return location.split('{')[0].strip()  # Remove any additional identifiers
    return location
# Function to process the metal-binding dataset
def process_metal_binding_data(file_path, output_path):
    # Read the dataset
    df = pd.read_csv(file_path)

    # Extract and clean the 'Subcellular_Location' column
    df['Subcellular_Location'] = df['Subcellular_Location'].apply(clean_location)

    # Extract the required columns
    df = df[['Metal_type', 'Subcellular_Location']]

    # Save the processed data with a more descriptive filename
    df.to_csv(os.path.join(output_path, 'metal_binding_subcellular_location.csv'), index=False)

# Function to process the non-metal binding dataset
def process_non_metal_binding_data(file_path, output_path):
    # Read the dataset
    df = pd.read_csv(file_path)

    # Extract and clean the 'Subcellular_Location' column
    df['Subcellular_Location'] = df['Subcellular_Location'].apply(clean_location)

    # We only need the 'Subcellular_Location' column with counts
    location_counts = df['Subcellular_Location'].value_counts().reset_index()
    location_counts.columns = ['Subcellular_Location', 'Count']

    # Save the processed data with a more descriptive filename
    location_counts.to_csv(os.path.join(output_path, 'non_metal_binding_subcellular_location_counts.csv'), index=False)

# Paths to the input files and output directory
metal_file_path = '.../refined_integrated_dataset_metals_final.csv'
nonmetal_file_path = '.../refined_integrated_dataset_nonmetal.csv'
output_dir = '.../'

# Process both datasets
process_metal_binding_data(metal_file_path, output_dir)
process_non_metal_binding_data(nonmetal_file_path, output_dir)
print("Data processing completed.")

#Great, let's plot this data and see where metal binding proteins are located:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap

# Function to read and process datasets
def read_and_process(file_path, metal_type=None):
    df = pd.read_csv(file_path)
    if metal_type:
        df['Metal_type'] = metal_type
    # Improve readability of 'Subcellular_Location' by replacing semicolons and commas
    df['Subcellular_Location'] = df['Subcellular_Location'].str.replace(';', ' |', regex=False)
    df['Subcellular_Location'] = df['Subcellular_Location'].str.replace(',', ' -', regex=False)
    return df[['Metal_type', 'Subcellular_Location']]

# Read the datasets
metal_df = read_and_process('.../metal_binding_subcellular_location.csv')
nonmetal_df = read_and_process('.../non_metal_binding_subcellular_location_counts.csv', 'non-metal')

# Combine the datasets
combined_df = pd.concat([metal_df, nonmetal_df], ignore_index=True)

# Aggregate the data by 'Metal_type' and 'Subcellular_Location' and count the occurrences
aggregated_df = combined_df.groupby(['Metal_type', 'Subcellular_Location']).size().reset_index(name='Count')

# Exclude the metal column called 'BR'
aggregated_df = aggregated_df[aggregated_df['Metal_type'] != 'BR']

def wrap_labels(labels, width):
    return ['\n'.join(textwrap.wrap(label, width)) for label in labels]

# Create a pivot table for the heatmap
pivot_df = aggregated_df.pivot(index='Subcellular_Location', columns='Metal_type', values='Count').fillna(0)

# We will use the log10 of the pivot_df for the actual heatmap
pivot_df_log = np.log10(pivot_df + 1)

# Plotting with black lines between cells
plt.figure(figsize=(18, 8))  # Adjust the figure size as necessary
ax = sns.heatmap(pivot_df_log, annot=False, cmap="YlGnBu", linecolor='black', linewidths=.5)

# Adjusting visual elements for readability
plt.xticks(rotation=45, fontsize=14, fontweight='bold')  # X-axis labels

# Wrapping y-axis labels
y_labels = pivot_df.index.tolist()
wrapped_y_labels = wrap_labels(y_labels, width=36)  # Adjust the width as necessary
ax.set_yticklabels(wrapped_y_labels, rotation=0, fontsize=11, fontweight='bold')

plt.xlabel('Metal Type', fontsize=14, fontweight='bold')  # X-axis title
plt.ylabel('Subcellular Location', fontsize=14, fontweight='bold')  # Y-axis title
# plt.title('Log-Scale Heatmap of Subcellular Location Abundance by Metal Type', fontsize=14, fontweight='bold')

# Update color bar labels with scientific notation
color_bar = ax.collections[0].colorbar
r = color_bar.vmax - color_bar.vmin
color_bar.set_ticks([color_bar.vmin + r / 4 * i for i in range(5)])
color_bar.set_ticklabels([f'$10^{{{int(i)}}}$' for i in range(5)])

# Add a title to the color bar to indicate the data is in log-scale
color_bar.ax.set_title('Count (log-scale)', pad=10, fontweight='bold')

plt.tight_layout()

# Save the figure to the desired path
plt.savefig('.../D_Subcellular_Location_Enhanced.png', dpi=300, bbox_inches='tight')
plt.show()
