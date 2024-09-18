import pandas as pd

# Load the dataset
df = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\enriched_binding_and_annotations_5.csv')

# Split the 'GO_Descriptions' into individual terms and explode the DataFrame
go_descriptions_series = df['GO_Descriptions'].str.split('; ').explode()

# Count the occurrences of each GO Term
go_descriptions_counts = go_descriptions_series.value_counts().nlargest(20)  # Top 10 GO Descriptions

import matplotlib.pyplot as plt
import seaborn as sns

def wrap_labels(labels, width):
    return ['\n'.join(textwrap.wrap(label, width)) for label in labels]

y_labels = go_descriptions_counts.index.tolist()
wrapped_y_labels = wrap_labels(y_labels, width=40)  # Adjust the width as necessary

# Plot for Top 10 GO Descriptions
plt.figure(figsize=(11, 8))
sns.barplot(x=go_descriptions_counts.values, y=wrapped_y_labels, color='purple')
# plt.title('Top 20 GO Term Descriptions', fontsize=14, fontweight='bold')
plt.xlabel('Frequency', fontsize=12, fontweight='bold')
plt.ylabel('GO Term Descriptions', fontsize=12, fontweight='bold')
plt.xticks(fontsize=8, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')

plt.tight_layout()
plt.show()
# Save the figure
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp2\\A_top_10_go_term_descriptions.png')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Split 'GO_Descriptions' into individual terms, explode, and then join back to include 'close_mutations_count'
go_terms_expanded = df.assign(GO_Descriptions=df['GO_Descriptions'].str.split('; ')).explode('GO_Descriptions')

# Calculate the mean 'close_mutations_count' for each GO term
mean_mutations_by_go = go_terms_expanded.groupby('GO_Descriptions')['close_mutations_count'].mean().reset_index()

# Sort by 'close_mutations_count' and get the top 10
top_mean_mutations_by_go = mean_mutations_by_go.nlargest(20, 'close_mutations_count')

# Plot
plt.figure(figsize=(13, 10))
barplot = sns.barplot(
    data=top_mean_mutations_by_go,
    y='GO_Descriptions',
    x='close_mutations_count',
    color='purple'
)

# Add value labels
for p in barplot.patches:
    barplot.annotate(
        '{:.2f}'.format(p.get_width()),
        (p.get_width(), p.get_y() + p.get_height() / 2),
        ha='left',
        va='center',
        fontsize=10,
        color='black',
        xytext=(5, 0),
        textcoords='offset points',
        fontweight='bold'
    )

#plt.title('Top 10 GO Term Descriptions by Mean Close Mutations Count')
plt.xlabel('Mean Close Mutations Count', fontsize=12, fontweight='bold')
plt.ylabel('GO Term Descriptions', fontsize=12, fontweight='bold')
plt.xticks(fontsize=14, fontweight='bold')
plt.yticks(fontsize=14, fontweight='bold')
plt.tight_layout()

plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp2\\B_top_10_go_term_descriptions_with_mean_mutation_count.png')
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
df = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\enriched_binding_and_annotations_5.csv')

# Group by 'Metal_type' and calculate the mean 'close_mutations_count'
metal_mutation_summary = df.groupby('Metal_type')['close_mutations_count'].mean().reset_index()

# Sort the summary DataFrame based on 'close_mutations_count'
metal_mutation_summary_sorted = metal_mutation_summary.sort_values('close_mutations_count', ascending=False)

# Create a bar plot for 'close_mutations_count' by 'Metal_type'
plt.figure(figsize=(14, 8))  # You can adjust the figure size as needed
barplot = sns.barplot(
    data=metal_mutation_summary_sorted,
    x='Metal_type',
    y='close_mutations_count',
    color='purple'
)
# Add value labels to each bar in the bar plot
for p in barplot.patches:
    barplot.annotate(
        format(p.get_height(), '.2f'),  # format the count with 2 decimal places
        (p.get_x() + p.get_width() / 2, p.get_height()),  # position for the text
        ha='center',  # center the text horizontally
        va='bottom',  # position the text at the bottom of the top side
        fontsize=10,  # font size
        color='black',  # text color
    fontweight = 'bold'
    )

#plt.title('Mean Close Mutations Count by Metal Type')
plt.xlabel('Metal Type', fontsize=12, fontweight='bold')
plt.ylabel('Mean Close Mutations Count', fontsize=12, fontweight='bold')
plt.xticks(fontsize=12, fontweight='bold', rotation=45)
plt.yticks(fontsize=10, fontweight='bold')
# Rotate the x-axis 12 for better readability
plt.tight_layout()  # Adjust the layout

# Save the figure
plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp2\\C_metal_type_vs_mutation_count.png')

plt.show()

# Fill missing Metal_type values with 'None' and ensure close_mutations_count is numeric
df['Metal_type'] = df['Metal_type'].fillna('None')
df['close_mutations_count'] = pd.to_numeric(df['close_mutations_count'], errors='coerce')

# Initialize binary columns for each metal to indicate its presence
all_metals = set(';'.join(df['Metal_type'].dropna().unique()).split(';'))
all_metals.discard('None')  # Remove 'None' from the set of metals

for metal in all_metals:
    df[metal] = df['Metal_type'].apply(lambda x: 1 if metal in x.split(';') else 0)

# Melt the DataFrame to have metal types as a single column
df_melted = pd.melt(df, id_vars=['CaseStudy_Category', 'close_mutations_count'],
                    value_vars=list(all_metals), var_name='Metal', value_name='Presence')

# Filter out rows where metal is not present
df_melted = df_melted[df_melted['Presence'] == 1]

# Visualization
plt.figure(figsize=(14, 8))
# Use the 'tab20' color palette for distinct colors
sns.swarmplot(data=df_melted, y='CaseStudy_Category', x='close_mutations_count',
              hue='Metal', size=14, palette='tab20', alpha=1)

#plt.title('Close Mutations Count by Case Study Category Highlighting Specific Metals', fontsize=14, fontweight='bold')
plt.ylabel('Case Study Category', fontsize=14, fontweight='bold')
plt.xlabel('Close Mutations Count', fontsize=12, fontweight='bold')
plt.xticks(fontsize=14, fontweight='bold')  # Rotate labels and align right
plt.yticks(fontsize=16, fontweight='bold')
plt.legend(title='Metal Type', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()

# Specify the path where you want to save the figure
figure_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp2\\D_overview_plot_with_distinct_colors.png'
plt.savefig(figure_path)  # Save the figure
print(f"Figure saved to {figure_path}")

plt.show()  # Display the figure

import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
df = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Updated_Targets_Casestudy_Unique.csv')

# Convert 'close_mutations_count' to a numeric data type, errors='coerce' will convert any problematic values to NaN
df['close_mutations_count'] = pd.to_numeric(df['close_mutations_count'], errors='coerce')

# Filter out rows with NaN in 'CaseStudy_Category' or 'close_mutations_count'
df_filtered = df.dropna(subset=['CaseStudy_Category', 'close_mutations_count'])

categories = df_filtered['CaseStudy_Category'].unique()
if len(categories) > 1:
    anova_data = [df_filtered[df_filtered['CaseStudy_Category'] == cat]['close_mutations_count'] for cat in categories]
    anova_result = stats.f_oneway(*anova_data)
    print(f"ANOVA result: F-statistic = {anova_result.statistic}, P-value = {anova_result.pvalue}")

    # Plotting the close_mutations_count for each CaseStudy_Category
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df_filtered, x='CaseStudy_Category', y='close_mutations_count', color='purple')
    #plt.title('Close Mutations Count across Case Study Categories')
    plt.ylabel('Close Mutations Count', fontsize=10, fontweight='bold')
    plt.xlabel('Case Study Category', fontsize=10, fontweight='bold')
    plt.xticks(rotation=45, fontsize=12, fontweight='bold') # Rotate the x-axis labels for better readability
    plt.yticks(fontsize=10, fontweight='bold')
    plt.tight_layout()  # Adjust the layout
    plt.savefig('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Supp2\\E_close_mutations_count_per_category.png')
    plt.show()
else:
    print("Not enough categories for ANOVA.")
