#################################
#######CASE STUDY 2: Stress#######
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
