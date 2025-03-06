#Figure 2
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pickle

# Load data
binding_df = pd.read_csv("binding_data.csv", low_memory=False)
mutations_per_metal_file = "mutations_per_metal.pkl"
pickle_file_path = "Distance_Mutations_MetalBinding.pkl"
mutations_near_metal_file = "Mutations_Near_Metal_Counts.pkl"

# First subplot: Histogram of distances (now wide across the figure)
def plot_histogram_of_distances(ax):
    binding_df_filtered = binding_df[binding_df['Distance'] <= 5]
    sns.histplot(binding_df_filtered['Distance'], bins=50, kde=True, color='skyblue', edgecolor='black', ax=ax)
    mean_distance = binding_df_filtered['Distance'].mean()
    median_distance = binding_df_filtered['Distance'].median()
    ax.axvline(mean_distance, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_distance:.2f} Å')
    ax.axvline(median_distance, color='green', linestyle='-', linewidth=2, label=f'Median: {median_distance:.2f} Å')
    ax.axvline(2, color='purple', linestyle='--', linewidth=2, label='Direct Coord (2 Å)')
    ax.axvline(3, color='purple', linestyle='--', linewidth=2, label='Coord Threshold (3 Å)')
    ax.axvline(5, color='orange', linestyle='--', linewidth=2, label='Cutoff (5 Å)')
    ax.set_xlim(0, 5)
    ax.set_ylabel('Frequency', fontsize=14, fontweight='bold')
    ax.set_xlabel('Euclidean Distance (Å) of Metal-Residue Interactions', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12, prop={'weight': 'bold'})  # Updated to use prop for fontweight

# Second subplot: Clustered heatmap for distances by metal type
def plot_clustered_heatmap(ax):
    thresholds = [0, 2, 3, 5]
    threshold_labels = ['<=2 Å', '2-3 Å', '3-5 Å']
    binding_df['Distance_Category'] = pd.cut(binding_df['Distance'], bins=thresholds, labels=threshold_labels, include_lowest=True)
    distance_counts = binding_df.groupby(['Metal_type', 'Distance_Category']).size().reset_index(name='Count')
    pivot_table = distance_counts.pivot(index="Distance_Category", columns="Metal_type", values="Count").fillna(0)
    pivot_table.replace(0, 0.1, inplace=True)  # Replace zeros to handle log scaling
    sns.heatmap(pivot_table, norm=LogNorm(), cmap="YlGnBu", ax=ax)
    ax.set_xlabel('Metal Type', fontsize=14, fontweight='bold')
    ax.set_ylabel('Distance Category', fontsize=14, fontweight='bold')

# Third subplot: Violin plot of distances by residue (with thinner bars)
def plot_violin_plot_of_distances(ax):
    amino_acids_3_letter = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    binding_df_residues = binding_df[binding_df['Residue_name'].isin(amino_acids_3_letter)].copy()  # Avoid SettingWithCopyWarning
    binding_df_residues.sort_values('Residue_name', inplace=True)
    sns.violinplot(data=binding_df_residues, x='Residue_name', y='Distance', hue='Residue_name', palette="husl", ax=ax, linewidth=1.0, legend=False)  # Set legend=False
    ax.set_xlabel('Residue Name', fontsize=14, fontweight='bold')
    ax.set_ylabel('Distance (Å)', fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=10)

# Fourth subplot: Histogram of mutation distances (swapped with the fifth)
def plot_histogram_of_mutation_distances(ax):
    with open(pickle_file_path, 'rb') as f:
        distances = pickle.load(f)
    ax.hist(distances, bins=30, color='purple', edgecolor='black')
    ax.set_xlabel('Distance to Metal-binding Site (Å)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Number of Mutations', fontsize=14, fontweight='bold')

# Fifth subplot: Bar plot of mutations near metal-binding sites (corrected)
def plot_metal_mutations_near_sites(ax):
    with open(mutations_near_metal_file, 'rb') as f:
        data = pickle.load(f)
    ale_count = data['ALE']
    ltee_count = data['LTEE']
    labels = ['ALE', 'LTEE']
    counts = [ale_count, ltee_count]
    ax.bar(labels, counts, color=['blue', 'green'])
    ax.set_xlabel('Dataset', fontsize=14, fontweight='bold')
    ax.set_ylabel('Mutation Count Near Metal (>1Å)', fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

# Create master figure layout
fig, axs = plt.subplots(3, 2, figsize=(12, 16), gridspec_kw={'height_ratios': [1, 1, 1]})
plt.subplots_adjust(hspace=0.2, wspace=0.3)  # Adjusted to reduce white space between subplots

# Adjust position of the first subplot to span across the top row
plot_histogram_of_distances(axs[0, 0])  # First row: wide subplot
axs[0, 0].set_position([0.05, 0.75, 0.9, 0.25])  # Adjusted to make it landscape across the top
axs[0, 1].axis('off')  # Empty subplot to fill space

plot_clustered_heatmap(axs[1, 0])  # Second row, left subplot
plot_violin_plot_of_distances(axs[1, 1])  # Second row, right subplot

plot_histogram_of_mutation_distances(axs[2, 0])  # Third row, left subplot
plot_metal_mutations_near_sites(axs[2, 1])  # Third row, right subplot (corrected)

# Adjust y-axis label positions to prevent overlap
axs[2, 1].yaxis.set_label_coords(1.1, 0.5)  # Move the y-axis label to the right

# Save the master figure
plt.savefig("Figure2\\Figure2.png", bbox_inches='tight')
plt.show()
