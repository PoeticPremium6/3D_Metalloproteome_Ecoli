#################################
#######CASE STUDY 2: Stress#######
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
df = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\enriched_binding_and_annotations_5.csv')

# Filter for Fe-binding proteins
df_fe = df[df['Metal_type'].str.contains('FE', na=False)]

# Split the 'GO_Descriptions' into individual terms and explode the DataFrame for frequency analysis
go_descriptions_series_fe = df_fe['GO_Descriptions'].str.split('; ').explode()
go_descriptions_counts_fe = go_descriptions_series_fe.value_counts().nlargest(10)

# Print the top 10 GO term descriptions by frequency
print("Top 10 GO Term Descriptions by Frequency:")
for term, count in go_descriptions_counts_fe.items():
    print(f"{term}: {count}")

# Prepare data for mean mutations count analysis
go_terms_expanded_fe = df_fe.assign(GO_Descriptions=df_fe['GO_Descriptions'].str.split('; ')).explode('GO_Descriptions')
mean_mutations_by_go_fe = go_terms_expanded_fe.groupby('GO_Descriptions')['close_mutations_count'].mean().reset_index()
top_mean_mutations_by_go_fe = mean_mutations_by_go_fe.nlargest(10, 'close_mutations_count')

# Print the top 10 GO term descriptions by mean close mutations count
print("\nTop 10 GO Term Descriptions by Mean Close Mutations Count:")
for index, row in top_mean_mutations_by_go_fe.iterrows():
    print(f"{row['GO_Descriptions']}: {row['close_mutations_count']}")

# Prepare data for mean mutations count analysis
go_terms_expanded_fe = df_fe.assign(GO_Descriptions=df_fe['GO_Descriptions'].str.split('; ')).explode('GO_Descriptions')
mean_mutations_by_go_fe = go_terms_expanded_fe.groupby('GO_Descriptions')['close_mutations_count'].mean().reset_index()
top_mean_mutations_by_go_fe = mean_mutations_by_go_fe.nlargest(10, 'close_mutations_count')

# Function to wrap labels
def wrap_labels(labels, width):
    return ['\n'.join(textwrap.wrap(label, width)) for label in labels]
# Create a figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(20, 10))

# Plot for Top 10 GO Descriptions by frequency
sns.barplot(ax=axes[0], x=go_descriptions_counts_fe.values, y=wrap_labels(go_descriptions_counts_fe.index, 20), color='purple')
#axes[0].set_title('Top 10 GO Term Descriptions by Frequency', fontsize=16, fontweight='bold')
axes[0].set_xlabel('Frequency', fontsize=14, fontweight='bold')
axes[0].set_ylabel('GO Term Descriptions', fontsize=14, fontweight='bold')
axes[0].tick_params(axis='x', labelsize=12,  labelrotation=0, width=2)
axes[0].tick_params(axis='y', labelsize=16, labelrotation=0, width=2)
plt.setp(axes[0].get_xticklabels(), fontweight='bold')
plt.setp(axes[0].get_yticklabels(), fontweight='bold')

# Plot for Top 10 GO Term Descriptions by Mean Close Mutations Count
sns.barplot(ax=axes[1], data=top_mean_mutations_by_go_fe, y=wrap_labels(top_mean_mutations_by_go_fe['GO_Descriptions'], 20), x='close_mutations_count', color='purple')
#axes[1].set_title('Top 10 GO Term Descriptions by Mean Close Mutations Count', fontsize=16, fontweight='bold')
axes[1].set_xlabel('Mean Close Mutations Count', fontsize=14, fontweight='bold')
axes[1].set_ylabel('', fontsize=14, fontweight='bold')  # Keep ylabel empty but adjust font size and weight for consistency
axes[1].tick_params(axis='x', labelsize=12, labelrotation=0, width=2)
axes[1].tick_params(axis='y', labelsize=16, labelrotation=0, width=2)
plt.setp(axes[1].get_xticklabels(), fontweight='bold')
plt.setp(axes[1].get_yticklabels(), fontweight='bold')

# Adjust layout
plt.tight_layout()
# Save the combined figure for Iron (Fe) analysis
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure5\\A_top_10_go_terms_combined_fe.png')
###

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the enriched dataset
df = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\enriched_binding_and_annotations_5.csv')

# Filter for Zinc-binding proteins
fe_df = df[df['Metal_type'] == 'FE']
hg_df = df[df['Metal_type'] == 'HG']

# Display the total number of Zinc-binding proteins identified
print(f"Total Fe-binding proteins identified: {fe_df.shape[0]}")
print(f"Total Hg-binding proteins identified: {hg_df.shape[0]}")

# Analyze mutation counts
# This step helps identify proteins with a significant number of mutations close to the metal-binding site, suggesting potential targets for further analysis or strain design.
fe_mutation_counts = fe_df.sort_values(by='close_mutations_count', ascending=False)
hg_mutation_counts = hg_df.sort_values(by='close_mutations_count', ascending=False)

# Display top entries for an overview
print(fe_mutation_counts[['UniProt_ID', 'Blattner Numbers', 'close_mutations_count']].head(10))
print(hg_mutation_counts[['UniProt_ID', 'Blattner Numbers', 'close_mutations_count']].head(10))

plt.figure(figsize=(10, 6))
sns.histplot(fe_df['close_mutations_count'], bins=30, kde=True, color='black')
plt.title('Distribution of Close Mutation Counts in Iron-Binding Proteins')
plt.xlabel('Close Mutation Count')
plt.ylabel('Frequency')
plt.grid(True)
plt.tight_layout()
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure5\\B_Distribution_of_Close_Mutation_Counts_Iron.png')
plt.show()

print(magnesium_mutation_counts[['UniProt_ID', 'Blattner Numbers', 'close_mutations_count']].head(10))
plt.figure(figsize=(10, 6))
sns.histplot(hg_df['close_mutations_count'], bins=30, kde=True, color='black')
plt.title('Distribution of Close Mutation Counts in Mercury-Binding Proteins')
plt.xlabel('Close Mutation Count')
plt.ylabel('Frequency')
plt.grid(True)
plt.tight_layout()
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp4\\B_Distribution_of_Close_Mutation_Counts_Mercury.png')
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the enriched dataset
df = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\enriched_binding_and_annotations_5.csv')

# Load the case study category dataset
case_study_df = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Updated_Targets_Casestudy_Unique.csv')

# Merge the enriched dataset with case study categories on Blattner Numbers
merged_df = pd.merge(df, case_study_df[['Blattner Numbers', 'CaseStudy_Category']], on='Blattner Numbers', how='left')

# Filter for Zinc (ZN) binding proteins
fe_df = merged_df[merged_df['Metal_type'].str.contains("FE", case=False, na=False)]
hg_df = merged_df[merged_df['Metal_type'].str.contains("HG", case=False, na=False)]

# Aggregate close mutation counts and count of zinc-binding sites for each Blattner Number
fe_agg = fe_df.groupby(['Blattner Numbers', 'CaseStudy_Category']).agg(
    metal_count=pd.NamedAgg(column='Metal_type', aggfunc='count'),
    total_close_mutations=pd.NamedAgg(column='close_mutations_count', aggfunc='sum')
).reset_index()
hg_agg = hg_df.groupby(['Blattner Numbers', 'CaseStudy_Category']).agg(
    metal_count=pd.NamedAgg(column='Metal_type', aggfunc='count'),
    total_close_mutations=pd.NamedAgg(column='close_mutations_count', aggfunc='sum')
).reset_index()

# Sort by metal_count for plotting
fe_agg_sorted =fe_agg.sort_values(by='metal_count', ascending=False)
hg_agg_sorted = hg_agg.sort_values(by='metal_count', ascending=False)

print("Analytics for Iron-binding proteins:")
print(fe_agg_sorted.head(10))

# Print analytics for Magnesium-binding proteins
print("\nAnalytics for Mercury-binding proteins:")
print(hg_agg_sorted.head(10))

# Plotting
plt.figure(figsize=(12, 8))
barplot = sns.barplot(x='Blattner Numbers', y='metal_count', hue='CaseStudy_Category', data=fe_agg_sorted.head(10), palette="viridis")

plt.yscale('log')  # Set the y-axis to logarithmic scale

# Add mutation count as annotation
for p in barplot.patches:
    height = p.get_height()
    plt.text(p.get_x() + p.get_width() / 2., height, f"{height:.1f}", ha='center', va='bottom')

plt.title('Top 10 Iron-Binding Proteins by Site Count with Mutation Counts, Colored by Case Study Category')
plt.xlabel('Protein/Gene (Blattner Numbers)')
plt.ylabel('Count of Iron-Binding Sites (Log Scale)')
plt.xticks(rotation=45)
plt.legend(title='Case Study Category')
plt.tight_layout()

# Save the figure - adjust the path as necessary
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure5\\C_Distribution_of_Iron_Binding_Sites_With_Mutation_Counts_Colored_By_CaseStudy_Category_Log_Scale.png')
plt.show()

# Plotting
plt.figure(figsize=(12, 8))
barplot = sns.barplot(x='Blattner Numbers', y='metal_count', hue='CaseStudy_Category', data=hg_agg_sorted.head(10), palette="viridis")

plt.yscale('log')  # Set the y-axis to logarithmic scale

# Add mutation count as annotation
for p in barplot.patches:
    height = p.get_height()
    plt.text(p.get_x() + p.get_width() / 2., height, f"{height:.1f}", ha='center', va='bottom')

plt.title('Top 10 Mercury-Binding Proteins by Site Count with Mutation Counts, Colored by Case Study Category')
plt.xlabel('Protein/Gene (Blattner Numbers)')
plt.ylabel('Count of Magnesium-Binding Sites (Log Scale)')
plt.xticks(rotation=45)
plt.legend(title='Case Study Category')
plt.tight_layout()

# Save the figure - adjust the path as necessary
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp4\\C_Distribution_of_Mercury_Binding_Sites_With_Mutation_Counts_Colored_By_CaseStudy_Category_Log_Scale.png')
plt.show()

####
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming 'data_path' and 'enriched_binding_data_path' are already defined and point to your CSV files
data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Updated_Targets_Casestudy_Unique.csv'
enriched_binding_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\enriched_binding_and_annotations_5.csv'

# Load the datasets
targets_df = pd.read_csv(data_path)
enriched_df = pd.read_csv(enriched_binding_data_path)

# Merge the datasets on 'Blattner Numbers'
merged_df = pd.merge(targets_df, enriched_df, on='Blattner Numbers', suffixes=('_targets', '_enriched'))

# Filter for Zinc (Zn) binding proteins in the merged dataset
fe_binding_proteins = merged_df[merged_df['Metal_type_enriched'].str.contains('FE', na=False)]
hg_binding_proteins = merged_df[merged_df['Metal_type_enriched'].str.contains('HG', na=False)]

# Split the 'ALE_Experiment' column into separate rows for each unique experiment
fe_binding_proteins_exploded = fe_binding_proteins.assign(ALE_Experiment=fe_binding_proteins['ALE_Experiment'].str.split(';')).explode('ALE_Experiment')
hg_binding_proteins_exploded = hg_binding_proteins.assign(ALE_Experiment=hg_binding_proteins['ALE_Experiment'].str.split(';')).explode('ALE_Experiment')

# Group by ALE experiments to see which ones are associated with high mutation counts in Zn-binding proteins
ale_impact = fe_binding_proteins_exploded.groupby('ALE_Experiment')['close_mutations_count_enriched'].sum().reset_index()
ale_impact = hg_binding_proteins_exploded.groupby('ALE_Experiment')['close_mutations_count_enriched'].sum().reset_index()

# Sort the results to highlight the ALE experiments with the highest total mutation counts in Zn-binding proteins
ale_impact_sorted = ale_impact.sort_values(by='close_mutations_count_enriched', ascending=False)

# Visualizing the impact of ALE experiments on Zn-binding proteins
plt.figure(figsize=(12, 8))
sns.barplot(x='ALE_Experiment', y='close_mutations_count_enriched', data=ale_impact_sorted, color='purple')
plt.title('Impact of ALE Experiments on Mutation Counts in Hg-Binding Proteins')
plt.xlabel('ALE Experiment')
plt.ylabel('Total Close Mutations Count')
plt.xticks(rotation=45)
plt.tight_layout()

# Save the figure
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp4\\D_Impact_of_ALE_Experiments_on_Hg_Proteins.png')
plt.show()
