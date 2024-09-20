import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import os

# Load the data
file_path = '.../PPI_Metal_Binding_Analysis.csv"
data = pd.read_csv(file_path)

# Drop rows with missing values in relevant columns
data.dropna(subset=['Degree', 'Mutation_Count', 'Log_Scaled_Mutations', 'Is_Metal_Binding'], inplace=True)

# Descriptive statistics for each group
stats_summary = data.groupby('Is_Metal_Binding').agg({
    'Degree': ['mean', 'std', 'count'],
    'Mutation_Count': ['mean', 'std'],
    'Log_Scaled_Mutations': ['mean', 'std']
})
print(stats_summary)

# Correlation analysis
metal_binding = data[data['Is_Metal_Binding'] == True]
non_metal_binding = data[data['Is_Metal_Binding'] == False]

correlation_metal = metal_binding[['Degree', 'Mutation_Count']].corr().iloc[0, 1]
correlation_non_metal = non_metal_binding[['Degree', 'Mutation_Count']].corr().iloc[0, 1]

print(f'Correlation (Metal Binding): {correlation_metal}')
print(f'Correlation (Non-Metal Binding): {correlation_non_metal}')

# Statistical tests
t_stat, p_value_ttest = stats.ttest_ind(
    metal_binding['Mutation_Count'].dropna(), 
    non_metal_binding['Mutation_Count'].dropna(), 
    equal_var=False
)

u_stat, p_value_mannwhitney = stats.mannwhitneyu(
    metal_binding['Mutation_Count'].dropna(), 
    non_metal_binding['Mutation_Count'].dropna()
)

print(f'T-test: T-statistic = {t_stat}, P-value = {p_value_ttest}')
print(f'Mann-Whitney U Test: U-statistic = {u_stat}, P-value = {p_value_mannwhitney}')

# Create output directory if it doesn't exist
output_dir = ".../"
os.makedirs(output_dir, exist_ok=True)

# Boxplot for mutation counts
plt.figure(figsize=(8, 6))
sns.boxplot(x='Is_Metal_Binding', y='Mutation_Count', data=data)
plt.title('Mutation Counts by Metal Binding Status')
plt.xlabel('Is Metal Binding')
plt.ylabel('Mutation Count')
plt.savefig(os.path.join(output_dir, 'mutation_counts_boxplot.png'))
plt.close()

# Scatter plot for Degree vs Mutation Count
plt.figure(figsize=(8, 6))
sns.scatterplot(x='Degree', y='Mutation_Count', hue='Is_Metal_Binding', data=data)
plt.title('Degree vs Mutation Count by Metal Binding Status')
plt.xlabel('Degree')
plt.ylabel('Mutation Count')
plt.legend(title='Is Metal Binding')
plt.savefig(os.path.join(output_dir, 'degree_vs_mutation_count_scatter.png'))
plt.close()

print("Figures saved successfully!")
