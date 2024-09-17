#################################
#######CASE STUDY 1: BASIC#######
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the enriched dataset
df = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\enriched_binding_and_annotations_5.csv')

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
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure4\\B_Distribution_of_Close_Mutation_Counts_Zinc.png')
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
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp3\\B_Distribution_of_Close_Mutation_Counts_Magnesium.png')
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
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure4\\C_Distribution_of_Zinc_Binding_Sites_With_Mutation_Counts_Colored_By_CaseStudy_Category_Log_Scale.png')
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
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp3\\C_Distribution_of_Magnesium_Binding_Sites_With_Mutation_Counts_Colored_By_CaseStudy_Category_Log_Scale.png')
plt.show()
