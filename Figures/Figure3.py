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
