#################################
#######CASE STUDY 1: BASIC#######
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap


# Load the dataset
df = pd.read_csv('.../enriched_binding_and_annotations_5.csv')

# Filter for Zn-binding proteins
df_zn = df[df['Metal_type'].str.contains('ZN', na=False)]

# Split the 'GO_Descriptions' into individual terms and explode the DataFrame for frequency analysis
go_descriptions_series_zn = df_zn['GO_Descriptions'].str.split('; ').explode()
go_descriptions_counts_zn = go_descriptions_series_zn.value_counts().nlargest(10)

# Print the top 10 GO term descriptions by frequency
print("Top 10 GO Term Descriptions by Frequency:")
for term, count in go_descriptions_counts_zn.items():
    print(f"{term}: {count}")

# Prepare data for mean mutations count analysis
go_terms_expanded_zn = df_zn.assign(GO_Descriptions=df_zn['GO_Descriptions'].str.split('; ')).explode('GO_Descriptions')
mean_mutations_by_go_zn = go_terms_expanded_zn.groupby('GO_Descriptions')['close_mutations_count'].mean().reset_index()
top_mean_mutations_by_go_zn = mean_mutations_by_go_zn.nlargest(10, 'close_mutations_count')

# Print the top 10 GO term descriptions by mean close mutations count
print("\nTop 10 GO Term Descriptions by Mean Close Mutations Count:")
for index, row in top_mean_mutations_by_go_zn.iterrows():
    print(f"{row['GO_Descriptions']}: {row['close_mutations_count']}")

# Prepare data for mean mutations count analysis
go_terms_expanded_zn = df_zn.assign(GO_Descriptions=df_zn['GO_Descriptions'].str.split('; ')).explode('GO_Descriptions')
mean_mutations_by_go_zn = go_terms_expanded_zn.groupby('GO_Descriptions')['close_mutations_count'].mean().reset_index()
top_mean_mutations_by_go_zn = mean_mutations_by_go_zn.nlargest(10, 'close_mutations_count')

# Function to wrap labels
def wrap_labels(labels, width):
    return ['\n'.join(textwrap.wrap(label, width)) for label in labels]
# Create a figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(20, 10))

# Plot for Top 10 GO Descriptions by frequency
sns.barplot(ax=axes[0], x=go_descriptions_counts_zn.values, y=wrap_labels(go_descriptions_counts_zn.index, 20), color='purple')
#axes[0].set_title('Top 10 GO Term Descriptions by Frequency', fontsize=16, fontweight='bold')
axes[0].set_xlabel('Frequency', fontsize=14, fontweight='bold')
axes[0].set_ylabel('GO Term Descriptions', fontsize=14, fontweight='bold')
axes[0].tick_params(axis='x', labelsize=12,  labelrotation=0, width=2)
axes[0].tick_params(axis='y', labelsize=16, labelrotation=0, width=2)
plt.setp(axes[0].get_xticklabels(), fontweight='bold')
plt.setp(axes[0].get_yticklabels(), fontweight='bold')

# Plot for Top 10 GO Term Descriptions by Mean Close Mutations Count
sns.barplot(ax=axes[1], data=top_mean_mutations_by_go_zn, y=wrap_labels(top_mean_mutations_by_go_zn['GO_Descriptions'], 20), x='close_mutations_count', color='purple')
#axes[1].set_title('Top 10 GO Term Descriptions by Mean Close Mutations Count', fontsize=16, fontweight='bold')
axes[1].set_xlabel('Mean Close Mutations Count', fontsize=14, fontweight='bold')
axes[1].set_ylabel('', fontsize=14, fontweight='bold')  # Keep ylabel empty but adjust font size and weight for consistency
axes[1].tick_params(axis='x', labelsize=12, labelrotation=0, width=2)
axes[1].tick_params(axis='y', labelsize=16, labelrotation=0, width=2)
plt.setp(axes[1].get_xticklabels(), fontweight='bold')
plt.setp(axes[1].get_yticklabels(), fontweight='bold')

# Adjust layout
plt.tight_layout()

# Save the combined figure
plt.savefig('.../A_top_10_go_terms_combined_zn.png')
plt.show()
###
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the enriched dataset
df = pd.read_csv('.../enriched_binding_and_annotations_5.csv')

# Filter for Zinc-binding proteins
zinc_df = df[df['Metal_type'] == 'ZN']
magnesium_df = df[df['Metal_type'] == 'MG']

# Display the total number of Zinc-binding proteins identified
print(f"Total Zinc-binding proteins identified: {zinc_df.shape[0]}")
print(f"Total Mg-binding proteins identified: {magnesium_df.shape[0]}")

# Analyze mutation counts
# This step helps identify proteins with a significant number of mutations close to the metal-binding site, suggesting potential targets for further analysis or strain design.
zinc_mutation_counts = zinc_df.sort_values(by='close_mutations_count', ascending=False)
magnesium_mutation_counts = magnesium_df.sort_values(by='close_mutations_count', ascending=False)

# Display top entries for Zinc-binding proteins
print("Top 10 Zinc-binding proteins by close mutations count:")
print(zinc_mutation_counts[['UniProt_ID', 'Blattner Numbers', 'close_mutations_count']].head(10))
# Display top entries for Magnesium-binding proteins
print("Top 10 Magnesium-binding proteins by close mutations count:")
print(magnesium_mutation_counts[['UniProt_ID', 'Blattner Numbers', 'close_mutations_count']].head(10))

plt.figure(figsize=(10, 6))
sns.histplot(zinc_df['close_mutations_count'], bins=30, kde=True, color='black')
#plt.title('Distribution of Close Mutation Counts in Zinc-Binding Proteins')
plt.xlabel('Close Mutation Count', fontsize=12, fontweight='bold')
plt.ylabel('Frequency', fontsize=12, fontweight='bold')
plt.xticks(fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
plt.grid(True)
plt.tight_layout()
plt.savefig('.../B_Distribution_of_Close_Mutation_Counts_Zinc.png')
plt.show()

print(magnesium_mutation_counts[['UniProt_ID', 'Blattner Numbers', 'close_mutations_count']].head(10))
plt.figure(figsize=(10, 6))
sns.histplot(magnesium_df['close_mutations_count'], bins=30, kde=True, color='black')
#plt.title('Distribution of Close Mutation Counts in Magnesium-Binding Proteins')
plt.xlabel('Close Mutation Count', fontsize=12, fontweight='bold')
plt.ylabel('Frequency', fontsize=12, fontweight='bold')
plt.xticks(fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
plt.grid(True)
plt.tight_layout()
plt.savefig('.../B_Distribution_of_Close_Mutation_Counts_Magnesium.png')
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the enriched dataset
df = pd.read_csv('.../enriched_binding_and_annotations_5.csv')

# Load the case study category dataset
case_study_df = pd.read_csv('.../Updated_Targets_Casestudy_Unique.csv')

# Merge the enriched dataset with case study categories on Blattner Numbers
merged_df = pd.merge(df, case_study_df[['Blattner Numbers', 'CaseStudy_Category']], on='Blattner Numbers', how='left')

# Filter for Zinc (ZN) binding proteins
zinc_df = merged_df[merged_df['Metal_type'].str.contains("ZN", case=False, na=False)]
magnesium_df = merged_df[merged_df['Metal_type'].str.contains("MG", case=False, na=False)]

# Aggregate close mutation counts and count of zinc-binding sites for each Blattner Number
zinc_agg = zinc_df.groupby(['Blattner Numbers', 'CaseStudy_Category']).agg(
    metal_count=pd.NamedAgg(column='Metal_type', aggfunc='count'),
    total_close_mutations=pd.NamedAgg(column='close_mutations_count', aggfunc='sum')
).reset_index()
magnesium_agg = magnesium_df.groupby(['Blattner Numbers', 'CaseStudy_Category']).agg(
    metal_count=pd.NamedAgg(column='Metal_type', aggfunc='count'),
    total_close_mutations=pd.NamedAgg(column='close_mutations_count', aggfunc='sum')
).reset_index()

# Sort by metal_count for plotting
zinc_agg_sorted = zinc_agg.sort_values(by='metal_count', ascending=False)
magnesium_agg_sorted = magnesium_agg.sort_values(by='metal_count', ascending=False)

# Print analytics for Zinc-binding proteins
print("Analytics for Zinc-binding proteins:")
print(zinc_agg_sorted.head(10))

# Print analytics for Magnesium-binding proteins
print("\nAnalytics for Magnesium-binding proteins:")
print(magnesium_agg_sorted.head(10))

# Plotting
plt.figure(figsize=(10, 8))  # Increase the width of the figure

# Create the bar plot
barplot = sns.barplot(
    x='Blattner Numbers',
    y='metal_count',
    hue='CaseStudy_Category',
    data=zinc_agg_sorted.head(10),
    palette="viridis",
    dodge=False)

plt.yscale('log')  # Set the y-axis to logarithmic scale
# Add mutation count as annotation with increased font size
for p in barplot.patches:
    height = p.get_height()
    barplot.annotate(f"{height:.1f}",
                     (p.get_x() + p.get_width() / 2., height),
                     ha='center', va='bottom',
                     fontsize=12, fontweight='bold')  # Adjust fontsize and fontweight

#plt.title('Top 10 Zinc-Binding Proteins by Site Count with Mutation Counts, Colored by Case Study Category', fontsize=16, fontweight='bold')
plt.xlabel('Protein/Gene (Blattner Numbers)', fontsize=14, fontweight='bold')
plt.ylabel('Count of Zinc-Binding Sites (Log Scale)', fontsize=14, fontweight='bold')
# Adjust x and y tick parameters for font size and boldness
plt.xticks(rotation=45, fontsize=14, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
# Adjust legend font size and title font size
plt.legend(title='Case Study Category', title_fontsize='13', fontsize='12', frameon=True, loc='upper right')

plt.tight_layout(pad=1.0)
plt.savefig('.../C_Distribution_of_Zinc_Binding_Sites_With_Mutation_Counts_Colored_By_CaseStudy_Category_Log_Scale.png')
plt.show()


# Plotting
plt.figure(figsize=(12, 8))  # Increase the width of the figure

# Create the bar plot
barplot = sns.barplot(
    x='Blattner Numbers',
    y='metal_count',
    hue='CaseStudy_Category',
    data=magnesium_agg_sorted.head(10),
    palette="viridis",
    dodge=False
)
plt.yscale('log')  # Set the y-axis to logarithmic scale

# Add mutation count as annotation with increased font size
for p in barplot.patches:
    height = p.get_height()
    barplot.annotate(f"{height:.1f}",
                     (p.get_x() + p.get_width() / 2., height),
                     ha='center', va='bottom',
                     fontsize=12, fontweight='bold')  # Adjust fontsize and fontweight

#plt.title('Top 10 Zinc-Binding Proteins by Site Count with Mutation Counts, Colored by Case Study Category', fontsize=16, fontweight='bold')
plt.xlabel('Protein/Gene (Blattner Numbers)', fontsize=8, fontweight='bold')
plt.ylabel('Count of Magnesium-Binding Sites (Log Scale)', fontsize=8, fontweight='bold')
# Adjust x and y tick parameters for font size and boldness
plt.xticks(rotation=45, fontsize=14, fontweight='bold', ha='right')  # Rotate labels and align right
plt.yticks(fontsize=10, fontweight='bold')
# Adjust legend font size and title font size
plt.legend(title='Case Study Category', title_fontsize='13', fontsize='12', frameon=True, loc='upper right')
# Save the figure - adjust the path as necessary
plt.savefig('.../C_Distribution_of_Magnesium_Binding_Sites_With_Mutation_Counts_Colored_By_CaseStudy_Category_Log_Scale.png')
plt.show()

# Assuming 'data_path' and 'enriched_binding_data_path' are already defined and point to your CSV files
data_path = 'C.../Updated_Targets_Casestudy_Unique.csv'
enriched_binding_data_path = '.../enriched_binding_and_annotations_5.csv'

# Load the datasets
targets_df = pd.read_csv(data_path)
enriched_df = pd.read_csv(enriched_binding_data_path)

# Merge the datasets on 'Blattner Numbers'
merged_df = pd.merge(targets_df, enriched_df, on='Blattner Numbers', suffixes=('_targets', '_enriched'))

# Filter for Zinc (Zn) binding proteins in the merged dataset
zn_binding_proteins = merged_df[merged_df['Metal_type_enriched'].str.contains('ZN', na=False)]
mg_binding_proteins = merged_df[merged_df['Metal_type_enriched'].str.contains('MG', na=False)]

# Split the 'ALE_Experiment' column into separate rows for each unique experiment
zn_binding_proteins_exploded = zn_binding_proteins.assign(ALE_Experiment=zn_binding_proteins['ALE_Experiment'].str.split(';')).explode('ALE_Experiment')
mg_binding_proteins_exploded = mg_binding_proteins.assign(ALE_Experiment=mg_binding_proteins['ALE_Experiment'].str.split(';')).explode('ALE_Experiment')

# Group by ALE experiments to see which ones are associated with high mutation counts in Zn-binding proteins
ale_impact = zn_binding_proteins_exploded.groupby('ALE_Experiment')['close_mutations_count_enriched'].sum().reset_index()
ale_impact = mg_binding_proteins_exploded.groupby('ALE_Experiment')['close_mutations_count_enriched'].sum().reset_index()

# Sort the results to highlight the ALE experiments with the highest total mutation counts in Zn-binding proteins
ale_impact_sorted = ale_impact.sort_values(by='close_mutations_count_enriched', ascending=False)

# Visualizing the impact of ALE experiments on Zn-binding proteins
plt.figure(figsize=(12, 8))
# Create barplot with purple bars
sns.barplot(x='ALE_Experiment', y='close_mutations_count_enriched', data=ale_impact_sorted, color='purple')
# Customize labels, ticks, and figure settings
plt.xlabel('ALE Experiment', fontsize=12, fontweight='bold')
plt.ylabel('Total Close Mutations Count', fontsize=12, fontweight='bold')
plt.xticks(rotation=45, fontsize=12, fontweight='bold', ha='right')
plt.yticks(fontsize=12, fontweight='bold')
plt.tight_layout()

# Save the figure
plt.savefig('.../D_Impact_of_ALE_Experiments_on_Zn_Proteins.png')
plt.show()

# Group by ALE experiments to see which ones are associated with high mutation counts in Zn-binding proteins
zn_ale_impact = zn_binding_proteins_exploded.groupby('ALE_Experiment')['close_mutations_count_enriched'].sum().reset_index()
zn_ale_impact_sorted = zn_ale_impact.sort_values(by='close_mutations_count_enriched', ascending=False)

# Group by ALE experiments to see which ones are associated with high mutation counts in Mg-binding proteins
mg_ale_impact = mg_binding_proteins_exploded.groupby('ALE_Experiment')['close_mutations_count_enriched'].sum().reset_index()
mg_ale_impact_sorted = mg_ale_impact.sort_values(by='close_mutations_count_enriched', ascending=False)

# Print analytics for Zn-binding proteins
print("Analytics for Zn-binding proteins:")
print(zn_ale_impact_sorted.head(10))

# Print analytics for Mg-binding proteins
print("\nAnalytics for Mg-binding proteins:")
print(mg_ale_impact_sorted.head(10))
