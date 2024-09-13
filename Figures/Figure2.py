#Let's move to Figure 2 which focuses on metals, mutations, residues, and proximity
#The plot displays the distances between metal ions and nearby residues in PDB structures
#X-Axis: Distance between metals and residues in Ångstroms.
#Y-Axis: Frequency of occurrence for each distance.
#Histogram Bars: Number of occurrences in each distance range.
#KDE Line: Smoothed representation of distance distribution.
#Key Insights:
#Peak Distances: Common interaction distances.
#Spread: Variation in interaction distances.
#Tail: Indicates residues farther from metals but still nearby.
#Average Distance: Central tendency of interactions.
#Outliers: Unusually large distances warranting investigation.
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the binding data from your local system
binding_df = pd.read_csv("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv")  # Update the path as needed

# Filter the dataframe for distances less than or equal to 5 Angstroms
binding_df_filtered = binding_df[binding_df['Distance'] <= 5]

# Count the frequency of interactions before and after 3 Å
total_interactions_before_3A = binding_df_filtered[binding_df_filtered['Distance'] <= 3].shape[0]
total_interactions_after_3A = binding_df_filtered[binding_df_filtered['Distance'] > 3].shape[0]

# Print the total frequency of interactions before and after 3 Å
print("Total interactions before 3 Å:", total_interactions_before_3A)
print("Total interactions after 3 Å:", total_interactions_after_3A)

# Set style
sns.set_style("whitegrid")

# Create the figure
plt.figure(figsize=(14, 8))

# Create the histogram with adjusted bins and a kernel density estimate
# This time, considering only data up to 5 Angstroms, spread equally
sns.histplot(binding_df_filtered['Distance'], bins=50, kde=True, color='skyblue', edgecolor='black')

# Calculate mean and median distances and add them as vertical lines for the filtered dataset
mean_distance = binding_df_filtered['Distance'].mean()
median_distance = binding_df_filtered['Distance'].median()
plt.axvline(mean_distance, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_distance:.2f} Å')
plt.axvline(median_distance, color='green', linestyle='-', linewidth=2, label=f'Median: {median_distance:.2f} Å')

# Add threshold lines for significant interaction distances within 5 Å
plt.axvline(2, color='purple', linestyle='--', linewidth=2, label='Direct Coordination Threshold (2 Å)')
plt.axvline(3, color='purple', linestyle='--', linewidth=2, label='Possible Coordination Threshold (3 Å)')
plt.axvline(5, color='orange', linestyle='--', linewidth=2, label='Interaction Cutoff (5 Å)')

# Add title and labels
#plt.title('Distribution of Distances (≤5 Å) Between Metals and Residues')
plt.ylabel('Frequency', fontsize=16, fontweight='bold')
plt.xlabel('Distance (Å)', fontsize=16, fontweight='bold')

# Adjusting the x-axis limit to make sure it's equally spread up to 5 Å
plt.xlim(0, 5)

# Customize x and y ticks
plt.xticks(fontsize=14, fontweight='bold')
plt.yticks(fontsize=14, fontweight='bold')

# Add a legend to the plot
legend = plt.legend(fontsize=14, title_fontsize=16)  # Add the legend
legend.get_title().set_fontweight('bold')  # Set the font weight for the legend title
# Set bold font properties for the legend text
for text in legend.get_texts():
    text.set_fontweight('bold')
    text.set_fontsize(14)
# Tighten the layout
plt.tight_layout()

# Save the figure
plt.savefig("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure2\\A_distance_distribution_equal_spread.png")  # Update the path as needed

# Show the plot
plt.show()

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd

# Load the binding data
binding_df = pd.read_csv("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv")  # Update the path as needed

# Define distance thresholds without the last one (5 Å)
thresholds = [0, 2, 3, 5]
threshold_labels = ['<=2 Å', '2-3 Å', '3-5 Å']

# Create a new column for distance categories based on thresholds
binding_df['Distance_Category'] = pd.cut(binding_df['Distance'],
                                         bins=thresholds,
                                         labels=threshold_labels,
                                         include_lowest=True)

# Aggregate data to count occurrences within each distance category for each metal type
distance_counts = binding_df.groupby(['Metal_type', 'Distance_Category']).size().reset_index(name='Count')

# Check the contents of distance_counts to ensure it was created correctly
print(distance_counts)

# Pivot the data to create a matrix for the heatmap
pivot_table = distance_counts.pivot(index="Distance_Category", columns="Metal_type", values="Count")

# Display the pivot table
print(pivot_table)

# Reorder columns based on hierarchical clustering distances
metal_types = pivot_table.columns.tolist()
clustered_metal_types = sns.clustermap(pivot_table[metal_types], method='complete', metric='euclidean', row_cluster=True, col_cluster=True).data2d.columns.tolist()
pivot_table = pivot_table.reindex(columns=clustered_metal_types)

# Flip the order of rows (Distance Category)
distance_categories = pivot_table.index.tolist()
pivot_table = pivot_table.reindex(distance_categories[::-1])

# Fill NaNs with 0
pivot_table = pivot_table.fillna(0)

# Replace zeros with a small positive value to handle log scaling
pivot_table.replace(0, 0.1, inplace=True)

# Determine the range of values in the pivot table
min_value = pivot_table.values.min()
max_value = pivot_table.values.max()

# Expand the color bar range to include another high value
vmin = min_value
vmax = max_value * 1.5  # Increase the maximum value by 20% to expand the color bar

# Create a clustermap with adjusted dendrogram ratio for x-axis
clustermap = sns.clustermap(
    pivot_table,
    cmap="YlGnBu",
    norm=LogNorm(vmin=vmin, vmax=vmax),  # Adjusted color bar range
    figsize=(18, 6),        # Larger figure size to accommodate dendrogram
    method='complete',      # Using 'complete' linkage method for better visibility of lower distance values
    metric='euclidean',     # Using 'euclidean' distance metric
    cbar_pos=(0.02, 0.85, 0.03, 0.1),  # Adjusting color bar position
    dendrogram_ratio=(0.1, 0.15),     # Adjusting dendrogram ratio, taller for y-axis, shorter for x-axis
    tree_kws={'linewidths': 1.5},    # Thicker lines for dendrogram
    row_cluster=True,               # Cluster rows based on distance
    col_cluster=True               # Cluster columns based on distance
)

# Adjust titles and labels with increased font size and boldness
clustermap.ax_heatmap.set_xlabel('Metal Type', fontsize=16, fontweight='bold', labelpad=10)  # Adjust label padding
clustermap.ax_heatmap.set_ylabel('Distance Category', fontsize=16, fontweight='bold', labelpad=10)  # Adjust label padding

# Calculate positions for centered ticks
xtick_positions = np.arange(pivot_table.shape[1])
ytick_positions = np.arange(pivot_table.shape[0])

# Ensure ticks and labels are properly shown
clustermap.ax_heatmap.set_xticks(xtick_positions)
clustermap.ax_heatmap.set_yticks(ytick_positions)

# Set tick labels explicitly and adjust alignment
clustermap.ax_heatmap.set_xticklabels(pivot_table.columns, rotation=45, ha='right', fontsize=13, fontweight='bold')
clustermap.ax_heatmap.set_yticklabels(pivot_table.index, fontsize=14, fontweight='bold', va='center')

# Add black grid lines between cells
for _, spine in clustermap.ax_heatmap.spines.items():
    spine.set_visible(True)  # Make all borders visible

clustermap.ax_heatmap.hlines(np.arange(1, pivot_table.shape[0]), *clustermap.ax_heatmap.get_xlim(), color='black', linewidth=1)
clustermap.ax_heatmap.vlines(np.arange(1, pivot_table.shape[1]), *clustermap.ax_heatmap.get_ylim(), color='black', linewidth=1)

# Remove major ticks
clustermap.ax_heatmap.tick_params(axis='both', which='major', length=0)

# Save the figure to the specified path with reduced whitespace
output_figure_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure2\\B_metal_distance_thresholds_clustered.png'
plt.savefig(output_figure_path, bbox_inches='tight', pad_inches=0.1)

#Violin Plot of Distances for Each Residue:
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the binding data
binding_df = pd.read_csv("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv")  # Update the path as needed

# Ensure only valid three-letter amino acid residues are considered
amino_acids_3_letter = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
binding_df = binding_df[binding_df['Residue_name'].isin(amino_acids_3_letter)]

# Sort the DataFrame by Residue_name to alphabetize the x-axis
binding_df.sort_values('Residue_name', inplace=True)

plt.figure(figsize=(14, 9))
# Use a vibrant and distinguishable palette, "husl"
palette = sns.color_palette("husl", len(amino_acids_3_letter))

# Create the violin plot with adjusted settings for broader violins
ax = sns.violinplot(data=binding_df, x='Residue_name', y='Distance', palette=palette,
                    linewidth=2.5, width=1.2, cut=0)

# Set the plot title and axis labels with improved font settings for readability and presentation
#plt.title('Alphabetized Distribution of Distances for Each Residue', fontsize=16, fontweight='bold')
plt.ylabel('Distance (Å)', fontsize=14, fontweight='bold')
plt.xlabel('Residue Name', fontsize=14, fontweight='bold')

# Add subtle vertical lines between categories for better visual separation
for x in np.arange(0.5, len(amino_acids_3_letter) - 0.5, 1):
    ax.axvline(x, color='gray', linestyle='--', alpha=0.6, linewidth=0.75)

# Rotate tick labels and set font size for better legibility
plt.xticks(rotation=45, fontsize=15, fontweight='bold')
plt.yticks(fontsize=15, fontweight='bold')

# Add a grid for easier data interpretation
plt.grid(True, linestyle='--', linewidth=0.6, alpha=0.5)

# Adjust plot layout to ensure complete label visibility
plt.subplots_adjust(bottom=0.15)

# Save the figure in high resolution for publication quality
plt.savefig("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure2\\C_residue_distance_violin.png", dpi=300, bbox_inches='tight')
plt.show()

#More residue plots
import pandas as pd

# Load the binding data
binding_df = pd.read_csv("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv")  # Update the path as needed

# Ensure only valid three-letter amino acid residues are considered
amino_acids_3_letter = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
binding_df = binding_df[binding_df['Residue_name'].isin(amino_acids_3_letter)]

# Group by metal type and residue name, then calculate statistical measures
aggregated_df = binding_df.groupby(['Metal_type', 'Residue_name'])['Distance'].agg(['min', 'max', 'mean', 'std']).reset_index()

# Rename columns for clarity
aggregated_df.columns = ['Metal_type', 'Residue_name', 'Min_Distance', 'Max_Distance', 'Mean_Distance', 'Std_Deviation']

# Save the aggregated data to a CSV file
csv_output_path = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\aggregated_metal_aa_interactions.csv"
aggregated_df.to_csv(csv_output_path, index=False)

print(f"Aggregated CSV file has been saved to: {csv_output_path}")

#Let's begin to study more about metal residue proximity, but integrate our mutation data
#Results will be used for Figures 2D/E and eventually Figure 3
import csv
import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import math

# Directories for input and output data
binding_data_file = "/mnt/data/project0032/Metalloproteome/Final/binding_data_New.csv"
mutation_files = [
    "/mnt/data/project0032/Metalloproteome/Final/Variants/ALE_Mutations.csv",
    "/mnt/data/project0032/Metalloproteome/Final/Variants/LTEE_Mutations.csv"
]
output_directory = "/mnt/data/project0032/Metalloproteome/Final/New_Output"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Load metal-binding data
def load_csv(file_path, delimiter=','):
    with open(file_path, 'r') as f:
        return [row for row in csv.DictReader(f, delimiter=delimiter)]

binding_data = load_csv(binding_data_file)

# Load mutation datasets (ALE and LTEE)
mutation_data = []
for mutation_file in mutation_files:
    mutation_data.extend(load_csv(mutation_file, delimiter=';'))

# Find matching mutations
matching_mutations = [
    {"mutation": mutation, "binding": binding}
    for mutation in mutation_data
    for binding in binding_data
    if mutation['AA_residue'] == binding['Residue_number']
]

# Count mutations per metal type
mutations_per_metal = defaultdict(int)
for entry in matching_mutations:
    mutations_per_metal[entry['binding']['Metal_type']] += 1

# Save the mutations per metal data to a pickle file
pickle_file = os.path.join(output_directory, "mutations_per_metal.pkl")
with open(pickle_file, 'wb') as f:
    pickle.dump(mutations_per_metal, f)

# Function to plot the data
def plot_data(df, save_path, title, xlabel, ylabel):
    plt.figure(figsize=(12, 8))
    plt.bar(df['Metal'], df['Count'], log=True)
    plt.xlabel(xlabel, fontsize=12, fontweight='bold')
    plt.ylabel(ylabel, fontsize=12, fontweight='bold')
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xticks(rotation=90, fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')
    plt.savefig(save_path)
    plt.show()

# Read and plot mutation per metal data
df = pd.DataFrame(list(mutations_per_metal.items()), columns=['Metal', 'Count']).sort_values(by='Count', ascending=False)
plot_save_path = os.path.join(output_directory, "Count_Metal_Mutations.png")
plot_data(df, plot_save_path, 'Descending Order Metal Counts on a Log Scale', 'Metal', 'Count (log scale)')

# Function to compute distance between coordinates
def compute_distance(coord):
    x, y, z = map(float, coord.replace("(", "").replace(")", "").split(","))
    return math.sqrt(x**2 + y**2 + z**2)

# Compute distances between mutations and metal-binding sites
distances = []
for mutation in mutation_data:
    closest_distance = float('inf')
    for binding in binding_data:
        if mutation['AA_residue'] == binding['Residue_number']:
            distance = compute_distance(binding['Residue_coord'])
            if distance < closest_distance:
                closest_distance = distance
    if closest_distance != float('inf'):
        distances.append(closest_distance)

# Save the distances to a pickle file
distance_pickle_file = os.path.join(output_directory, "Distance_Mutations_MetalBinding.pkl")
with open(distance_pickle_file, 'wb') as f:
    pickle.dump(distances, f)

# Plot the histogram of distances
plt.figure(figsize=(10, 6))
plt.hist(distances, bins=30, color='purple', edgecolor='black')
plt.xlabel('Distance from Reference Point to Metal-binding Site (Ångström)', fontsize=14, fontweight='bold')
plt.ylabel('Number of Mutations', fontsize=14, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.5)
plt.xticks(fontsize=14, fontweight='bold')
plt.yticks(fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(output_directory, "D_Distance_Mutations_MetalBinding.png"), dpi=300)
plt.show()

# Count mutations near metal-binding sites (distance threshold of 5Å)
def mutations_near_metal(metalloproteins, mutations, threshold=5):
    count = 0
    for metal in metalloproteins:
        metal_residue = int(metal['Residue_number'])
        for mutation in mutations:
            if abs(metal_residue - int(mutation['AA_residue'])) <= threshold:
                count += 1
    return count

ale_count = mutations_near_metal(binding_data, [m for m in mutation_data if 'ALE' in m['Dataset']])
ltee_count = mutations_near_metal(binding_data, [m for m in mutation_data if 'LTEE' in m['Dataset']])

# Save mutation counts near metal-binding sites
mutation_count_pickle = os.path.join(output_directory, "Mutations_Near_Metal_Counts.pkl")
with open(mutation_count_pickle, 'wb') as f:
    pickle.dump({'ALE': ale_count, 'LTEE': ltee_count}, f)

# Plot the counts
labels = ['ALE', 'LTEE']
counts = [ale_count, ltee_count]
plt.bar(labels, counts, color=['blue', 'green'])
plt.xlabel('Dataset', fontsize=14, fontweight='bold')
plt.ylabel('Mutation Count Near Metal (≤ 5Å)', fontsize=14, fontweight='bold')
plt.xticks(fontsize=14, fontweight='bold')
plt.yticks(fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(output_directory, "E_Mutations_Near_Metal.png"))
plt.show()

