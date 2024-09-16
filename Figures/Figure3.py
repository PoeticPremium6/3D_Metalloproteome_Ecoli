#5. Comparative Analysis of Mutation rate of proteins  Metal-Binding and Non-Metal Binding Proteins
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Specify the necessary columns for analysis
metal_columns = ['UniProt_ID', 'Gene_Length', 'Metal_type']  # Add any other relevant columns you need
non_metal_columns = ['UniProt_ID', 'Gene_Length']  # Adjust for non-metal data

# Load only necessary columns from the datasets
metal_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\refined_integrated_dataset_metals_final.csv'
non_metal_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\refined_integrated_dataset_nonmetal.csv'
mutation_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\categorized_dataset.csv'

metal_data = pd.read_csv(metal_data_path, usecols=metal_columns, low_memory=False)
non_metal_data = pd.read_csv(non_metal_data_path, usecols=non_metal_columns, low_memory=False)
mutation_data = pd.read_csv(mutation_data_path, low_memory=False)

# Rename 'Entry' to 'UniProt_ID' in mutation_data for consistency
mutation_data.rename(columns={'Entry': 'UniProt_ID'}, inplace=True)

# Count the number of mutations for each UniProt ID
mutation_counts = mutation_data['UniProt_ID'].value_counts().reset_index()
mutation_counts.columns = ['UniProt_ID', 'Number_of_Mutations']

# Merge the mutation data with metal and non-metal data
metal_mutations = pd.merge(metal_data, mutation_counts, on='UniProt_ID', how='left')
non_metal_mutations = pd.merge(non_metal_data, mutation_counts, on='UniProt_ID', how='left')

# Merge Grantham_Category from the mutation data
metal_mutations = pd.merge(metal_mutations, mutation_data[['UniProt_ID', 'Grantham_Category']], on='UniProt_ID', how='left')
non_metal_mutations = pd.merge(non_metal_mutations, mutation_data[['UniProt_ID', 'Grantham_Category']], on='UniProt_ID', how='left')

# Calculate mutation rates (mutations per gene length)
metal_mutations['Mutation_Rate'] = metal_mutations['Number_of_Mutations'] / metal_mutations['Gene_Length']
non_metal_mutations['Mutation_Rate'] = non_metal_mutations['Number_of_Mutations'] / non_metal_mutations['Gene_Length']

# Group by Grantham_Category and Metal_Type, then calculate average mutation rate
mutation_rate_metal = metal_mutations.groupby(['Grantham_Category', 'Metal_type'])['Mutation_Rate'].mean().unstack(fill_value=0)

# Add the non-metal data as a separate column
non_metal_avg_rate = non_metal_mutations.groupby('Grantham_Category')['Mutation_Rate'].mean()
mutation_rate_metal['Non-Metal'] = non_metal_avg_rate

# Apply log scale transformation
mutation_rate_metal_log = np.log1p(mutation_rate_metal)

# Flip the y-axis by sorting the index in descending order
mutation_rate_metal_log = mutation_rate_metal_log.sort_index(ascending=False)

# Create the heatmap with gridlines and log scale using the reversed 'viridis' colormap for correct color mapping
plt.figure(figsize=(12, 8))
ax = sns.heatmap(mutation_rate_metal_log, annot=True, cmap='viridis_r', fmt=".2f",
                 cbar_kws={'label': 'Log-scaled Average Mutation Rate'},
                 linewidths=0.5, linecolor='lightgrey')  # Add gridlines

plt.title('Heatmap of Log-scaled Average Mutation Rate by Grantham Category and Metal Type')
plt.xlabel('Metal Type', fontsize=14, fontweight='bold')
plt.ylabel('Grantham Category', fontsize=14, fontweight='bold')
plt.xticks(rotation=45)
plt.yticks(rotation=0)  # Ensure the y-axis labels are horizontal

# Save the figure (commented out for now)
 plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure3\\A_Mutation_Rate_Heatmap_Log_Scale.png')
# Show the heatmap
plt.tight_layout()
plt.show()

#Plot AA Property Change and Grantham Score (Figure 3B)
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\processed_merged_data_distance_5.csv'
output_dir = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure3\\"
distance = 5.0

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
def plot_mutations(data_path, output_dir, distance):
    # Load the dataset
    df = pd.read_csv(data_path)

    # Ensure 'Grantham_Score' is numeric
    df['Grantham_Score'] = pd.to_numeric(df['Grantham_Score'], errors='coerce')

    # Apply the categorization function to Grantham scores
    df['Grantham_Category'] = df['Grantham_Score'].apply(categorize_score)

    # Drop rows with missing values in these columns
    df.dropna(subset=['Property_Change', 'Grantham_Score', 'Metal_type'], inplace=True)

    # Sort the dataframe
    df.sort_values(by=['Metal_type', 'Property_Change'], inplace=True)

    # Define colors for each Grantham category
    grantham_colors = {
        'conservative': 'green',
        'moderately conservative': 'blue',
        'moderately radical': 'yellow',
        'radical': 'red',
        'unknown': 'grey'
    }

    # Create a plot resembling a heatmap
    plt.figure(figsize=(12, 8))
    sns.scatterplot(data=df, x='Metal_type', y='Property_Change', hue='Grantham_Category', palette=grantham_colors,
                    s=200)
    plt.title(f'Metal Type vs Property Change Categorized by Grantham Score (Distance ≤ {distance})')
    plt.xlabel('Metal Type')
    plt.ylabel('Property Change')
    plt.xticks(rotation=45)

    # Add legend
    plt.legend(title='Grantham Category', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Save the plot
    plt.savefig(f"{output_dir}metalloprotein_heatmap_style_plot_distance_{distance}.png", bbox_inches='tight')
    plt.close()
# Define the output directory

# File paths for different distances
file_paths = {
    5: "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\processed_merged_data_distance_5.csv",
    10: "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\processed_merged_data_distance_10.csv",
    20: "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\processed_merged_data_distance_20.csv"
}

# Generate plots for each file
for distance, file_path in file_paths.items():
    plot_mutations(file_path, output_dir, distance)

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
pickle_file_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\filtered_enhanced_data.pkl'
data = pd.read_pickle(pickle_file_path)

# Convert 'impact' to boolean for easier calculations
data['impact'] = data['impact'].astype(bool)

# Function to calculate proportion of disruption
def calculate_disruption_proportion(group):
    return group['impact'].mean()

# Analysis by Transcription Factor
disruption_by_tf = data.groupby('tf').apply(calculate_disruption_proportion)
disruption_by_tf = disruption_by_tf.reset_index(name='Disruption_Proportion')

# Analysis by Transcription Factor and Mutation Category (Experimental)
disruption_tf_experimental = data.groupby(['tf', 'Mutation_Category_Experimental']).apply(calculate_disruption_proportion)
disruption_tf_experimental = disruption_tf_experimental.reset_index(name='Disruption_Proportion')

# Pivot the data to create a matrix for the heatmap
heatmap_data = disruption_tf_experimental.pivot("Mutation_Category_Experimental", "tf", "Disruption_Proportion")
heatmap_data.fillna(0, inplace=True)  # Fill NaN values with zero

# Reverse the order of the y-axis
heatmap_data = heatmap_data.sort_index(ascending=False)

# Analysis by Metal Type
disruption_metal = data.groupby('Metal_type').apply(calculate_disruption_proportion)
disruption_metal = disruption_metal.reset_index(name='Disruption_Proportion')

# Define directory to save figures
figures_dir = '/Users/josspa/GPS-M/Figures/4.0/Figure3'

# Plotting
sns.set(style="whitegrid")

# Adjust the figure size and aspect ratio
fig, ax = plt.subplots(figsize=(28, 10))  # Increased width for better readability

# Plot for Disruption by TF and Mutation Category (Experimental) as Heatmap
heatmap = sns.heatmap(heatmap_data, ax=ax, cmap="viridis_r", cbar_kws={'label': 'Proportion of Disruption'}, linewidths=.5)
#plt.title('Disruption Proportion by TF and Experimental Mutation Category')
plt.xlabel('Transcription Factor', fontsize=14, fontweight='bold')
plt.ylabel('Mutation Category (Experimental)', fontsize=14, fontweight='bold')
plt.xticks(rotation=45, fontsize=16, fontweight='bold')
plt.yticks(rotation=0, fontsize=15, fontweight='bold')
plt.tight_layout()  # Adjust layout to fit everything properly
# Customizing colorbar (legend) label
cbar = heatmap.collections[0].colorbar
cbar.set_label('Proportion of Disruption', fontsize=14, fontweight='bold')

plt.savefig(f"{figures_dir}/C_Disruption_by_TF_Experimental_Mutation_Category_Heatmap.png")
plt.close()

# Plot for Metal Type
disruption_metal_sorted = disruption_metal.sort_values('Disruption_Proportion', ascending=False)

# Assuming 'disruption_metal_sorted' is a predefined list or array
figures_dir = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure3\\'

# Create barplot with single color (purple)
plt.figure(figsize=(12, 8))
bar = sns.barplot(x='Metal_type', y='Disruption_Proportion', data=disruption_metal_sorted, color='purple')
# Customize plot appearance
plt.xlabel('Metal Type', fontsize=15, fontweight='bold')
plt.ylabel('Proportion of Disruption', fontsize=16, fontweight='bold')
plt.xticks(rotation=45, fontsize=18, fontweight='bold')
plt.yticks(rotation=0, fontsize=14, fontweight='bold')
# Save the updated plot
plt.savefig(f"{figures_dir}/D_Disruption_Metal_Type_purple.png")
plt.close()

#4. Analysis of Post-Translational Modifications (PTMs) in the Context of Mutations
import pandas as pd
import matplotlib.pyplot as plt

# Load up datasets
base_directory = '/Users/josspa/GPS-M/'
binding_data = pd.read_csv(base_directory + 'binding_data_New.csv')
non_metal_data = pd.read_csv(base_directory + 'non_metal_data.csv')
ptm_data = pd.read_csv(base_directory + 'MutFunc/other_ptms.tab', sep='\t')  # Assuming it's tab-delimited
mutation_data = pd.read_csv(base_directory + 'categorized_dataset.csv')

# Merge PTM data with binding and mutation data on UniProt ID
binding_ptms = pd.merge(binding_data, ptm_data, left_on='UniProt_ID', right_on='acc', how='inner')
non_metal_ptms = pd.merge(non_metal_data, ptm_data, left_on='UniProt_ID', right_on='acc', how='inner')
mutation_ptms = pd.merge(mutation_data, ptm_data, left_on='Uniprot_ID', right_on='acc', how='inner')

# Count the PTMs for each dataset
binding_ptm_counts = binding_ptms['modification'].value_counts().reset_index(name='Binding')
non_metal_ptm_counts = non_metal_ptms['modification'].value_counts().reset_index(name='Non-Binding')
mutation_ptm_counts = mutation_ptms['modification'].value_counts().reset_index(name='Mutation')

# Merge the counts into a single DataFrame
merged_counts = binding_ptm_counts.merge(non_metal_ptm_counts, on='index', how='outer')
merged_counts = merged_counts.merge(mutation_ptm_counts, on='index', how='outer').fillna(0)

# Sort the data for visualization
merged_counts.sort_values(by='Binding', ascending=False, inplace=True)

# Plotting
fig, ax = plt.subplots(figsize=(14, 8))
# Width of the bars
#bar_width = 0.35

# Set position of bar on X axis
r1 = range(len(merged_counts))
r2 = [x + bar_width + 0.05 for x in r1]
r3 = [x + 2 * (bar_width + 0.05) for x in r1]

# Make the plot
ax.bar(r1, merged_counts['Binding'], color='blue', width=bar_width, edgecolor='black', label='Metal-Binding')
ax.bar(r2, merged_counts['Non-Binding'], color='green', width=bar_width, edgecolor='black', label='Non-Metal-Binding')
ax.bar(r3, merged_counts['Mutation'], color='orange', width=bar_width, edgecolor='black', label='Mutation')

# Add xticks on the middle of the group bars
ax.set_xlabel('Modification Type', fontsize=12, fontweight='bold')
# Set xticks and labels
bar_width = 0.30  # Adjust as per your bar width setup
ax.set_xticks([r + bar_width for r in range(len(merged_counts))])
# Boldface y-axis tick labels
ax.set_xticklabels(merged_counts['index'], rotation=90, fontsize=10, fontweight='bold')  # Increased font size for xtick labels
# Set ylabel with increased font size and bold font
ax.set_ylabel('Count', fontsize=12, fontweight='bold')
# Set yscale to logarithmic
ax.set_yscale('log')
# Set title with increased font size and bold font
#ax.set_title('Distribution of PTMs in Metal-Binding vs Non-Metal-Binding vs Mutation Data', fontsize=14, fontweight='bold')

# Create legend & Show graphic
ax.legend()

# Ensure layout is tight so labels are not cut off
plt.tight_layout()

# Save the plot
save_path = base_directory + 'Figures/4.0/Figure3/E_PTM_Distribution_Comparison_Grouped.png'
plt.savefig(save_path, format='png', dpi=300)
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
file_path = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Final\\Results\\PPI_Metal_Binding_Analysis.csv"
data = pd.read_csv(file_path)

# Drop duplicate rows
data = data.drop_duplicates()

# Replace missing values with 0
data['Mutation_Count'] = data['Mutation_Count'].fillna(0)
data['Log_Scaled_Mutations'] = data['Log_Scaled_Mutations'].fillna(0)

# Strip any extra spaces from column names
data.columns = data.columns.str.strip()

# Filter out rows where 'Mutation_Count' or 'Log_Scaled_Mutations' are zero
filtered_data = data[(data['Mutation_Count'] > 0) & (data['Log_Scaled_Mutations'] > 0)]

# Group by 'Binding_Metal_Types' and aggregate
grouped_data = filtered_data.groupby('Binding_Metal_Types').agg({
    'Mutation_Count': 'sum',
    'Log_Scaled_Mutations': 'mean',
    'Degree': ['sum', 'mean'],  # Sum and mean of degrees
    'UniProt_ID': 'count'  # Counting the number of proteins
}).rename(columns={'UniProt_ID': 'Protein_Count'}).reset_index()

# Flatten MultiIndex columns
grouped_data.columns = ['_'.join(col).strip() if isinstance(col, tuple) else col for col in grouped_data.columns.values]
grouped_data = grouped_data.rename(columns={
    'Degree_sum': 'Total_Degree',
    'Degree_mean': 'Average_Degree'
})

# Prepare the plot
plt.figure(figsize=(14, 8))

# Create a bar plot for Log-Scaled Mutations in purple
sns.barplot(x='Binding_Metal_Types_', y='Log_Scaled_Mutations_mean', data=grouped_data, color='purple', label='Avg Log-Scaled Mutations')

# Create a line plot for Total Degree
sns.lineplot(x='Binding_Metal_Types_', y='Total_Degree_sum', data=grouped_data, color='orange', marker='o', label='Total Degree')

# Add labels and title
plt.xlabel('Binding Metal Type')
plt.ylabel('Average Log-Scaled Mutations')
plt.title('Average Log-Scaled Mutations and Total Degree by Metal Binding Type')
plt.xticks(rotation=45)

# Create a secondary y-axis for Degree
ax2 = plt.gca().twinx()
ax2.set_ylabel('Total Degree', color='orange')
ax2.tick_params(axis='y', labelcolor='orange')

# Combine y-axis labels
ax1 = plt.gca()
ax1.set_ylabel('Average Log-Scaled Mutations', color='purple')
ax1.tick_params(axis='y', labelcolor='purple')

# Add legend
plt.legend(loc='upper left')

# Save the plot
plt.savefig("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure3\\Integrated_Plot.png", bbox_inches='tight', dpi=300)
plt.close()
