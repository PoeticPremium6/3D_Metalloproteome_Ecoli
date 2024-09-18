#################################
#######CASE STUDY 3: Industrial#######
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
df = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\enriched_binding_and_annotations_5.csv')

# Filter for MO-binding proteins
df_mo = df[df['Metal_type'].str.contains('MO', na=False)]

# Split the 'GO_Descriptions' into individual terms and explode the DataFrame for frequency analysis
go_descriptions_series_mo = df_mo['GO_Descriptions'].str.split('; ').explode()
go_descriptions_counts_mo = go_descriptions_series_mo.value_counts().nlargest(10)

# Print the top 10 GO term descriptions by frequency
print("Top Molybdenum 10 GO Term Descriptions by Frequency:")
for term, count in go_descriptions_counts_mo.items():
    print(f"{term}: {count}")
(f"{row['GO_Descriptions']}: {row['close_mutations_count']}")

# Prepare data for mean mutations count analysis
go_terms_expanded_mo = df_mo.assign(GO_Descriptions=df_mo['GO_Descriptions'].str.split('; ')).explode('GO_Descriptions')
mean_mutations_by_go_mo = go_terms_expanded_mo.groupby('GO_Descriptions')['close_mutations_count'].mean().reset_index()
top_mean_mutations_by_go_mo = mean_mutations_by_go_mo.nlargest(10, 'close_mutations_count')

# Function to wrap labels
def wrap_labels(labels, width):
    return ['\n'.join(textwrap.wrap(label, width)) for label in labels]
# Create a figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(20, 10))

# Plot for Top 10 GO Descriptions by frequency
sns.barplot(ax=axes[0], x=go_descriptions_counts_mo.values, y=wrap_labels(go_descriptions_counts_mo.index, 20), color='purple')
#axes[0].set_title('Top 10 GO Term Descriptions by Frequency', fontsize=16, fontweight='bold')
axes[0].set_xlabel('Frequency', fontsize=14, fontweight='bold')
axes[0].set_ylabel('GO Term Descriptions', fontsize=14, fontweight='bold')
axes[0].tick_params(axis='x', labelsize=12,  labelrotation=0, width=2)
axes[0].tick_params(axis='y', labelsize=16, labelrotation=0, width=2)
plt.setp(axes[0].get_xticklabels(), fontweight='bold')
plt.setp(axes[0].get_yticklabels(), fontweight='bold')

# Plot for Top 10 GO Term Descriptions by Mean Close Mutations Count
sns.barplot(ax=axes[1], data=top_mean_mutations_by_go_mo, y=wrap_labels(top_mean_mutations_by_go_mo['GO_Descriptions'], 20), x='close_mutations_count', color='purple')
#axes[1].set_title('Top 10 GO Term Descriptions by Mean Close Mutations Count', fontsize=16, fontweight='bold')
axes[1].set_xlabel('Mean Close Mutations Count', fontsize=14, fontweight='bold')
axes[1].set_ylabel('', fontsize=14, fontweight='bold')  # Keep ylabel empty but adjust font size and weight for consistency
axes[1].tick_params(axis='x', labelsize=12, labelrotation=0, width=2)
axes[1].tick_params(axis='y', labelsize=16, labelrotation=0, width=2)
plt.setp(axes[1].get_xticklabels(), fontweight='bold')
plt.setp(axes[1].get_yticklabels(), fontweight='bold')

# Adjust layout
plt.tight_layout()

# Save the combined figure for Molybdenum (MO) analysis
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure6\\A_top_10_go_terms_combined_mo.png')
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the enriched dataset
df = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\enriched_binding_and_annotations_5.csv')

# Filter for Zinc-binding proteins
mo_df = df[df['Metal_type'] == 'MO']
pb_df = df[df['Metal_type'] == 'PB']

# Display the total number of Zinc-binding proteins identified
print(f"Total Mo-binding proteins identified: {mo_df.shape[0]}")
print(f"Total Pb-binding proteins identified: {pb_df.shape[0]}")

# Analyze mutation counts
# This step helps identify proteins with a significant number of mutations close to the metal-binding site, suggesting potential targets for further analysis or strain design.
mo_mutation_counts = mo_df.sort_values(by='close_mutations_count', ascending=False)
pb_mutation_counts = pb_df.sort_values(by='close_mutations_count', ascending=False)

# Display top entries for Zinc-binding proteins
print("Top 10 Zinc-binding proteins by close mutations count:")
print(mo_mutation_counts[['UniProt_ID', 'Blattner Numbers', 'close_mutations_count']].head(10))
# Display top entries for Magnesium-binding proteins
print("Top 10 Magnesium-binding proteins by close mutations count:")
print(pb_mutation_counts[['UniProt_ID', 'Blattner Numbers', 'close_mutations_count']].head(10))

plt.figure(figsize=(10, 6))
sns.histplot(mo_df['close_mutations_count'], bins=30, kde=True, color='black')
#plt.title('Distribution of Close Mutation Counts in Zinc-Binding Proteins')
plt.xlabel('Close Mutation Count', fontsize=12, fontweight='bold')
plt.ylabel('Frequency', fontsize=12, fontweight='bold')
plt.xticks(fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
plt.grid(True)
plt.tight_layout()
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure6\\B_Distribution_of_Close_Mutation_Counts_Molybdenum.png')
plt.show()

plt.figure(figsize=(10, 6))
sns.histplot(pb_df['close_mutations_count'], bins=30, kde=True, color='black')
#plt.title('Distribution of Close Mutation Counts in Magnesium-Binding Proteins')
plt.xlabel('Close Mutation Count', fontsize=12, fontweight='bold')
plt.ylabel('Frequency', fontsize=12, fontweight='bold')
plt.xticks(fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
plt.grid(True)
plt.tight_layout()
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp3\\B_Distribution_of_Close_Mutation_Counts_Lead.png')

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
mo_df = merged_df[merged_df['Metal_type'].str.contains("MO", case=False, na=False)]
pb_df = merged_df[merged_df['Metal_type'].str.contains("PB", case=False, na=False)]

# Aggregate close mutation counts and count of zinc-binding sites for each Blattner Number
mo_agg = mo_df.groupby(['Blattner Numbers', 'CaseStudy_Category']).agg(
    metal_count=pd.NamedAgg(column='Metal_type', aggfunc='count'),
    total_close_mutations=pd.NamedAgg(column='close_mutations_count', aggfunc='sum')
).reset_index()
pb_agg = pb_df.groupby(['Blattner Numbers', 'CaseStudy_Category']).agg(
    metal_count=pd.NamedAgg(column='Metal_type', aggfunc='count'),
    total_close_mutations=pd.NamedAgg(column='close_mutations_count', aggfunc='sum')
).reset_index()

# Sort by metal_count for plotting
mo_agg_sorted = mo_agg.sort_values(by='metal_count', ascending=False)
pb_agg_sorted = pb_agg.sort_values(by='metal_count', ascending=False)

print("Analytics for Zinc-binding proteins:")
print(mo_agg_sorted.head(10))

# Print analytics for Magnesium-binding proteins
print("\nAnalytics for Magnesium-binding proteins:")
print(pb_agg_sorted.head(10))
# Plotting
plt.figure(figsize=(10, 8))  # Increase the width of the figure
# Create the bar plot
barplot = sns.barplot(
    x='Blattner Numbers',
    y='metal_count',
    hue='CaseStudy_Category',
    data=mo_agg_sorted.head(10),
    palette="viridis",
    dodge=False)

plt.yscale('log')  # Set the y-axis to logarithmic scale
# Add mutation count as annotation
for p in barplot.patches:
    height = p.get_height()
    barplot.annotate(f"{height:.1f}",
                     (p.get_x() + p.get_width() / 2., height),
                     ha='center', va='bottom',
                     fontsize=12, fontweight='bold')  # Adjust fontsize and fontweight

#plt.title('Top 10 Molybdenum-Binding Proteins by Site Count with Mutation Counts, Colored by Case Study Category', fontsize=16, fontweight='bold')
plt.xlabel('Protein/Gene (Blattner Numbers)', fontsize=14, fontweight='bold')
plt.ylabel('Count of Molybdenum-Binding Sites (Log Scale)', fontsize=14, fontweight='bold')
# Adjust x and y tick parameters for font size and boldness
plt.xticks(rotation=45, fontsize=14, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
# Adjust legend font size and title font size
plt.legend(title='Case Study Category', title_fontsize='13', fontsize='12', frameon=True, loc='upper right')
plt.tight_layout(pad=1.0)
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure6\\C_Distribution_of_Molybdenum_Binding_Sites_With_Mutation_Counts_Colored_By_CaseStudy_Category_Log_Scale.png')
plt.show()

# Plotting
plt.figure(figsize=(10, 8))  # Increase the width of the figure
# Create the bar plot
barplot = sns.barplot(
    x='Blattner Numbers',
    y='metal_count',
    hue='CaseStudy_Category',
    data=pb_agg_sorted.head(10),
    palette="viridis",
    dodge=False)

plt.yscale('log')  # Set the y-axis to logarithmic scale

# Add mutation count as annotation
for p in barplot.patches:
    height = p.get_height()
    barplot.annotate(f"{height:.1f}",
                     (p.get_x() + p.get_width() / 2., height),
                     ha='center', va='bottom',
                     fontsize=12, fontweight='bold')  # Adjust fontsize and fontweight

#plt.title('Top 10 Lead-Binding Proteins by Site Count with Mutation Counts, Colored by Case Study Category', fontsize=16, fontweight='bold')
plt.xlabel('Protein/Gene (Blattner Numbers)', fontsize=14, fontweight='bold')
plt.ylabel('Count of Lead-Binding Sites (Log Scale)', fontsize=14, fontweight='bold')
# Adjust x and y tick parameters for font size and boldness
plt.xticks(rotation=45, fontsize=14, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
# Adjust legend font size and title font size
plt.legend(title='Case Study Category', title_fontsize='13', fontsize='12', frameon=True, loc='upper right')

plt.tight_layout(pad=1.0)
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp5\\C_Distribution_of_Lead_Binding_Sites_With_Mutation_Counts_Colored_By_CaseStudy_Category_Log_Scale.png')
plt.show()

###
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
mo_binding_proteins = merged_df[merged_df['Metal_type_enriched'].str.contains('MO', na=False)]
pb_binding_proteins = merged_df[merged_df['Metal_type_enriched'].str.contains('PB', na=False)]

# Split the 'ALE_Experiment' column into separate rows for each unique experiment
mo_binding_proteins_exploded = mo_binding_proteins.assign(ALE_Experiment=mo_binding_proteins['ALE_Experiment'].str.split(';')).explode('ALE_Experiment')
pb_binding_proteins_exploded = pb_binding_proteins.assign(ALE_Experiment=pb_binding_proteins['ALE_Experiment'].str.split(';')).explode('ALE_Experiment')

# Group by ALE experiments to see which ones are associated with high mutation counts in Zn-binding proteins
ale_impact = mo_binding_proteins_exploded.groupby('ALE_Experiment')['close_mutations_count_enriched'].sum().reset_index()
ale_impact = pb_binding_proteins_exploded.groupby('ALE_Experiment')['close_mutations_count_enriched'].sum().reset_index()

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
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure6\\D_Impact_of_ALE_Experiments_on_Mo_Proteins.png')
plt.show()
