#Figure 3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy import stats

# Load data
input_csv_path = 'mutation_rate_metal_log.csv'
mutation_rate_metal_log = pd.read_csv(input_csv_path, index_col=0)
data_path = 'processed_merged_data_distance_5.csv'
df = pd.read_csv(data_path)
binding_data = pd.read_csv("binding_data.csv")
non_metal_data = pd.read_csv("non_metal_data.csv")
ptm_data = pd.read_csv("MutFunc\\other_ptms.tab", sep='\t')
data = pd.read_csv("CPPI_Metal_Binding_Analysis.csv")

# Output directory
output_dir = 'Figure3'
os.makedirs(output_dir, exist_ok=True)

# Preprocess correlation data
data.dropna(subset=['Degree', 'Mutation_Count', 'Is_Metal_Binding'], inplace=True)
metal_binding = data[data['Is_Metal_Binding'] == True]
non_metal_binding = data[data['Is_Metal_Binding'] == False]

# Calculate correlations
correlation_metal, p_value_metal = stats.pearsonr(metal_binding['Degree'], metal_binding['Mutation_Count'])
correlation_non_metal, p_value_non_metal = stats.pearsonr(non_metal_binding['Degree'], non_metal_binding['Mutation_Count'])

# Create figure with 3 rows and 2 columns
fig = plt.figure(figsize=(20, 18))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.2], hspace=0.4, wspace=0.4)

# ---- Subplot 1: Heatmap ----
ax1 = fig.add_subplot(gs[0, 0])
sns.heatmap(
    mutation_rate_metal_log,
    ax=ax1,
    annot=True,
    cmap='viridis_r',
    fmt=".2f",
    annot_kws={"size": 8, "weight": 'bold'},
    cbar_kws={'label': 'Log-scaled Average Mutation Rate'},
    linewidths=0.5,
    linecolor='lightgrey',
    xticklabels=True,
    yticklabels=True
)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right', fontsize=12, fontweight='bold')
ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0, fontsize=12, fontweight='bold')
#ax1.set_title('Heatmap of Log-scaled Average Mutation Rate', fontsize=16, fontweight='bold')
ax1.set_xlabel('Metal Type', fontsize=14, fontweight='bold')
ax1.set_ylabel('Grantham Category', fontsize=14, fontweight='bold')

# ---- Subplot 2: Scatter Plot ----
ax2 = fig.add_subplot(gs[0, 1])
grantham_colors = {
    'conservative': 'green',
    'moderately conservative': 'blue',
    'moderately radical': 'yellow',
    'radical': 'red',
    'unknown': 'grey'
}
sns.scatterplot(data=df, x='Metal_type', y='Property_Change', hue='Grantham_Category', palette=grantham_colors,
                s=200, ax=ax2)
#ax2.set_title('Metal Type vs Property Change (Distance ≤ 5.0)', fontsize=16, fontweight='bold')
ax2.set_xlabel('Metal Type', fontsize=14, fontweight='bold')
ax2.set_ylabel('Property Change', fontsize=14, fontweight='bold')

# ---- Subplot 3: Metal-Binding Correlation ----
ax3 = fig.add_subplot(gs[1, 0])
sns.regplot(x='Degree', y='Mutation_Count', data=metal_binding, ax=ax3,
            scatter_kws={'color': 'purple'}, line_kws={'color': 'purple'})
#ax3.set_title('Correlation: Metal-Binding Proteins', fontsize=16, fontweight='bold')
ax3.set_xlabel('Degree', fontsize=14, fontweight='bold')
ax3.set_ylabel('Mutation Count', fontsize=14, fontweight='bold')
ax3.text(0.6, 0.85, f'Correlation: {correlation_metal:.2f}\nP-value: {p_value_metal:.2e}',
         fontsize=12, color='purple', fontweight='bold', transform=ax3.transAxes)

# ---- Subplot 4: Non-Metal-Binding Correlation ----
ax4 = fig.add_subplot(gs[1, 1])
sns.regplot(x='Degree', y='Mutation_Count', data=non_metal_binding, ax=ax4,
            scatter_kws={'color': 'gray'}, line_kws={'color': 'gray'})
#ax4.set_title('Correlation: Non-Metal-Binding Proteins', fontsize=16, fontweight='bold')
ax4.set_xlabel('Degree', fontsize=14, fontweight='bold')
ax4.set_ylabel('Mutation Count', fontsize=14, fontweight='bold')
ax4.text(0.6, 0.85, f'Correlation: {correlation_non_metal:.2f}\nP-value: {p_value_non_metal:.2e}',
         fontsize=12, color='gray', fontweight='bold', transform=ax4.transAxes)

# ---- Subplot 5: PTM Modification (Landscape Format) ----
ax5 = fig.add_subplot(gs[2, :])  # Span entire row
common_ids_binding = set(binding_data['UniProt_ID']) & set(ptm_data['acc'])
common_ids_non_metal = set(non_metal_data['UniProt_ID']) & set(ptm_data['acc'])
common_rows_binding = binding_data[binding_data['UniProt_ID'].isin(common_ids_binding)]
common_rows_non_metal = non_metal_data[non_metal_data['UniProt_ID'].isin(common_ids_non_metal)]
merged_common_binding = pd.merge(common_rows_binding, ptm_data, left_on='UniProt_ID', right_on='acc', how='left')
merged_common_non_metal = pd.merge(common_rows_non_metal, ptm_data, left_on='UniProt_ID', right_on='acc', how='left')
mod_count_binding = merged_common_binding['modification'].value_counts().reset_index()
mod_count_binding.columns = ['modification', 'Binding']
mod_count_non_binding = merged_common_non_metal['modification'].value_counts().reset_index()
mod_count_non_binding.columns = ['modification', 'Non-Binding']
merged_counts = pd.merge(mod_count_binding, mod_count_non_binding, on='modification', how='outer').fillna(0)
merged_counts = merged_counts.sort_values(by='Binding', ascending=False)
bar_width = 0.4
x = range(len(merged_counts))
ax5.bar(x, merged_counts['Binding'], width=bar_width, color='purple', edgecolor='black', label='Binding')
ax5.bar([pos + bar_width for pos in x], merged_counts['Non-Binding'], width=bar_width, color='gray', edgecolor='black', label='Non-Binding')
ax5.set_xticks([pos + bar_width / 2 for pos in x])
ax5.set_xticklabels(merged_counts['modification'], rotation=45, ha="right", fontsize=12, fontweight='bold')
ax5.set_ylabel('Count (Log Scale)', fontsize=14, fontweight='bold')
ax5.set_xlabel('Modification Type', fontsize=14, fontweight='bold')
ax5.set_yscale('log')
ax5.legend(fontsize=12, title='Category', title_fontsize=12)
#ax5.set_title('PTM Modification Distribution', fontsize=16, fontweight='bold')

# Adjust layout
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'Figure_3.png'), bbox_inches='tight', dpi=300)
plt.show()
