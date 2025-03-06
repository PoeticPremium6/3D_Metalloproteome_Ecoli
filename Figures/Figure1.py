import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import numpy as np
import textwrap

# Load your main datasets
binding_df = pd.read_csv("binding_data.csv", low_memory=False)
non_metal_df = pd.read_csv("non_metal_data.csv")
metals_data_path = 'integrated_dataset_metals_final.csv'
non_metals_data_path = 'integrated_dataset_nonmetal.csv'
ec_class_path = 'EC_Class.csv'

# Load additional datasets
metals_data = pd.read_csv(metals_data_path, low_memory=False)
non_metals_data = pd.read_csv(non_metals_data_path, low_memory=False)

# Label the data
metals_data['Type'] = 'Metal'
non_metals_data['Type'] = 'Non-Metal'
non_metals_data['Metal_type'] = 'Non-Metal'  # Add 'Metal_type' column for consistency

# Combine the datasets
combined_data = pd.concat([metals_data, non_metals_data], ignore_index=True)

# Extracting the first part of the EC number and format it to match EC_Class
combined_data['EC_Number_Modified'] = combined_data['EC_Number'].apply(
    lambda x: str(x).split('.')[0] + '.-.-.-' if pd.notnull(x) else None)

# Determine the top-4 metals
top_metal_types = combined_data[combined_data['Type'] == 'Metal']['Metal_type'].value_counts().head()

# Filter the data to include only top-4 metals and non-metals
filtered_data = combined_data[
    combined_data['Metal_type'].isin(top_metal_types.index) | (combined_data['Type'] == 'Non-Metal')]

# Mapping modified EC numbers to their classes
ec_class_data = pd.read_csv(ec_class_path, low_memory=False)
mapped_data = pd.merge(filtered_data, ec_class_data, left_on='EC_Number_Modified', right_on='EC_Number', how='left')

# Data Analysis
grouped_data = mapped_data.groupby(['Metal_type', 'Class', 'Class_Function']).size().reset_index(name='Count')

# Create the figure with 3 rows and 2 columns
fig, axs = plt.subplots(3, 2, figsize=(25, 15))  # Increase figure size
plt.subplots_adjust(hspace=0.5, wspace=0.4)  # Adjusted spacing to reduce white space

# 1A: Pie chart of metal-binding vs metal-free proteins
labels = ['Metal-binding', 'Metal-free']
sizes = [len(binding_df['PDB_name'].unique()), len(non_metal_df['PDB_name'].unique())]
colors = ['#7f7fff', '#ff7f7f']

# Create the pie chart with bold and larger text
wedges, texts, autotexts = axs[0, 0].pie(
    sizes,
    labels=labels,
    colors=colors,
    autopct='%1.1f%%',
    startangle=90,
    textprops={'fontsize': 16, 'fontweight': 'bold'}  # Bold and larger text for labels
)

# Adjust the properties of the percentage text
for autotext in autotexts:
    autotext.set_fontsize(14)  # Set the font size for percentages
    autotext.set_fontweight('bold')  # Make the percentages bold

axs[0, 0].axis('equal')
# axs[0, 0].set_title("Metal-binding vs Metal-free Proteins", fontsize=20, fontweight='bold')

# 1B: Pie chart for the count of metalloproteins by metal type
metal_protein_counts = combined_data[combined_data['Type'] == 'Metal']['Metal_type'].value_counts()
top_10_metal_types = metal_protein_counts.head(10)

sizes = list(top_10_metal_types.values)
labels = list(top_10_metal_types.index)

colors = plt.cm.tab20(np.linspace(0, 1, len(labels)))

# Create the pie chart without percentages
wedges, texts = axs[0, 1].pie(sizes, colors=colors, startangle=90, textprops={'fontsize': 20, 'fontweight': 'bold'})
axs[0, 1].axis('equal')
#axs[0, 1].set_title("Top 10 Metalloproteins by Type", fontsize=24, fontweight='bold')

# Prepare legend data
percentages = [f'{size / sum(sizes) * 100:.1f}%' for size in sizes]  # Calculate percentages
legend_labels = [f'{label} ({percent})' for label, percent in zip(labels, percentages)]

# Add the legend to the right of the pie chart
axs[0, 1].legend(legend_labels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20, title='Metal Types')

plt.tight_layout()  # Adjust layout to fit larger pie chart

# 2A: Log-scaled distribution of metalloproteins
axs[1, 0].bar(metal_protein_counts.index, np.log10(metal_protein_counts.values + 1), color='purple')
#axs[1, 0].set_title("Log-Scaled Distribution of Metalloproteins", fontsize=20, fontweight='bold')
axs[1, 0].set_xlabel("Metal Type", fontsize=16, fontweight='bold')  # Made bold
axs[1, 0].set_ylabel("Number of Metalloproteins (Log-Scale)", fontsize=16, fontweight='bold')  # Made bold
axs[1, 0].tick_params(axis='x', rotation=45, labelsize=12)

# Adjust x and y axis ticks for the bar plot
for tick in axs[1, 0].get_xticklabels():
    tick.set_fontsize(14)
    tick.set_fontweight('bold')

for tick in axs[1, 0].get_yticklabels():
    tick.set_fontsize(14)
    tick.set_fontweight('bold')

# 2B: Placeholder for GO terms (empty)
#axs[1, 1].set_title("Placeholder for GO Terms", fontsize=20, fontweight='bold')
axs[1, 1].axis('off')  # Turn off axis for placeholder

# 3A: Heatmap of metabolic pathways
combined_data['Metabolic_Pathway_Short'] = combined_data['Metabolic_Pathway'].apply(
    lambda x: x.split(';')[0].split('.')[0] if isinstance(x, str) else None
)

pivot_df = pd.pivot_table(combined_data, index="Metabolic_Pathway_Short", columns="Metal_type", aggfunc="size", fill_value=0)

sns.heatmap(np.log10(pivot_df + 1), ax=axs[2, 0], cmap="YlGnBu", linewidths=1.5, linecolor='black', cbar=False)
#axs[2, 0].set_title("Log-scaled Heatmap of Metabolic Pathways", fontsize=20, fontweight='bold')
axs[2, 0].set_xlabel("Metal Type", fontsize=16, fontweight='bold')  # Made bold
axs[2, 0].set_ylabel("Metabolic Pathway", fontsize=18, fontweight='bold')  # Made bold

# Adjust x and y axis ticks for the metabolic pathways heatmap
for tick in axs[2, 0].get_xticklabels():
    tick.set_fontsize(14)
    tick.set_fontweight('bold')  # Make x-axis ticks bold

for tick in axs[2, 0].get_yticklabels():
    tick.set_fontsize(14)
    tick.set_fontweight('bold')

# 3B: Heatmap of protein families
# Read the protein family data
protein_family_df = pd.read_csv('protein_family_data.csv')

# Replace missing values
protein_family_df['Protein_Families'] = protein_family_df['Protein_Families'].fillna('NA')

# Aggregate and count occurrences
aggregated_df = protein_family_df.groupby(['Metal_type', 'Protein_Families']).size().reset_index(name='Count')

# Select the top 10 protein families
aggregated_total_counts = aggregated_df.groupby('Protein_Families')['Count'].sum().reset_index()
top_10_protein_families = aggregated_total_counts.sort_values('Count', ascending=False).head(10)['Protein_Families']

# Filter the original aggregated_df to only include these top 10 protein families
filtered_df = aggregated_df[aggregated_df['Protein_Families'].isin(top_10_protein_families)]

# Create pivot table for heatmap of the top 10 protein families
pivot_df_proteins = filtered_df.pivot(index='Protein_Families', columns='Metal_type', values='Count').fillna(0)
pivot_df_proteins_log = np.log10(pivot_df_proteins + 1)

# Plotting protein families heatmap
sns.heatmap(pivot_df_proteins_log, ax=axs[2, 1], annot=False, cmap="YlGnBu", linewidths=0.5, linecolor='black', cbar=False)

# Adjust title and labels
#axs[2, 1].set_title("Log-scaled Heatmap of Protein Families", fontsize=20, fontweight='bold')
axs[2, 1].set_xlabel("Metal Type", fontsize=16, fontweight='bold')  # Made bold
axs[2, 1].set_ylabel("Protein Families", fontsize=16, fontweight='bold')  # Made bold

# Adjust x and y axis ticks for the protein families heatmap
for tick in axs[2, 1].get_xticklabels():
    tick.set_fontsize(14)
    tick.set_fontweight('bold')  # Make x-axis ticks bold

axs[2, 1].set_yticklabels(pivot_df_proteins.index, rotation=0, horizontalalignment='left', fontsize=16, fontweight='bold')
axs[2, 1].yaxis.tick_right()  # Move y-axis ticks to the right side
axs[2, 1].yaxis.set_label_position("right")  # Move the label to the right

plt.tight_layout()

# Save the figure to the specified path
plt.savefig(
    'Figure1\\Master_Figure.png',
    dpi=300
)
plt.show()

# Let's plot a description of the top-5 GO term Descriptions for each category for metals vs. non-metals
# Function to get the top n GO terms for a given category and data
metals_data_path = 'merged_metals_data.csv'

# Read the data
metals_data = pd.read_csv(metals_data_path)

def get_top_n_go_terms(data, category, n=5):
    top_terms = (
        data[data['Category'] == category]
            .groupby('Description')
            .size()
            .reset_index(name='Count')
            .sort_values(by='Count', ascending=False)
            .head(n)
    )
    return top_terms

# Define the categories you want to plot
categories = ['molecular_function', 'biological_process', 'cellular_component']

# Function to wrap text
def wrap_labels(ax, width, axis='y'):
    if axis == 'y':
        labels = ax.get_yticklabels()
    elif axis == 'x':
        labels = ax.get_xticklabels()

    wrapped_labels = [textwrap.fill(label.get_text(), width) for label in labels]

    if axis == 'y':
        ax.set_yticklabels(wrapped_labels, rotation=0, ha='right')
    elif axis == 'x':
        ax.set_xticklabels(wrapped_labels, rotation=0, ha='right')

# Start plotting
fig, axes = plt.subplots(3, 1, figsize=(25, 15))  # Adjusted for one column for metal-binding proteins only
plt.subplots_adjust(hspace=0.5)  # Adjusted vertical spacing

for i, category in enumerate(categories):
    # Get the top 5 GO term descriptions for metal-binding genes
    top_metal_terms = get_top_n_go_terms(metals_data, category, n=5)

    # Plot for metal-binding genes (Purple color)
    ax1 = sns.barplot(y='Description', x='Count', data=top_metal_terms, ax=axes[i], color='purple')
    ax1.set_title(f'{category} - Metal-binding', fontweight='bold', fontsize=26)  # Increased Title Font Size
    ax1.set_xscale('log')
    ax1.set_xlabel('Count (Log Scale)', fontweight='bold', fontsize=22)  # Increased X-label Font Size
    ax1.set_ylabel('', fontweight='bold', fontsize=20)

    # Set x-axis tick labels properties
    ax1.set_xticks([1, 10, 100, 1000])
    ax1.get_xaxis().set_major_formatter(plt.ScalarFormatter(useMathText=True))
    for tick in ax1.get_xticklabels():
        tick.set_fontweight('bold')
        tick.set_fontsize(26)  # Increased X-tick Font Size

    # Set y-axis tick labels properties
    for tick in ax1.get_yticklabels():
        tick.set_fontweight('bold')
        tick.set_fontsize(32)  # Increased Y-tick Font Size

    # Set tick parameters for bold tick marks
    ax1.tick_params(axis='x', which='both', width=2.0)
    ax1.tick_params(axis='y', which='both', width=2.0)

    # Wrap y-tick labels
    wrap_labels(ax1, width=26, axis='y')

# Adjust layout with increased padding to prevent clipping of labels
plt.tight_layout(pad=5.0)  # Adjusted padding
plt.show()


# Save the figure
output_figure_path = 'Figure1\\SP_2B_Go_terms_metal_only.png'
plt.savefig(output_figure_path, dpi=300)  # Save at higher resolution
