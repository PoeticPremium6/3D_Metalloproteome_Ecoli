#Let's make a plot to count how many times metals interact with residues
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Path to your CSV file
csv_file_path = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv"

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
output_figure_path = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp1\\A_amino_acid_residue_distribution.png"
plt.savefig(output_figure_path, bbox_inches='tight', pad_inches=0.1)

print(f"Figure saved to: {output_figure_path}")

#Now let's make a plot to show the amount of times each metal interacts with a amino acid resiude
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy  # Ensure scipy is installed

# Path to your CSV file
csv_file_path = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv"

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
output_figure_path = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp1\\B_amino_acid_metal_interaction_counts_clustermap.png"
plt.savefig(output_figure_path, bbox_inches='tight', pad_inches=0.1)
