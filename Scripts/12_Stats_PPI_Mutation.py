import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import os

# Load the data
file_path = .../PPI_Metal_Binding_Analysis.csv"
data = pd.read_csv(file_path)

# Drop rows with missing values
data.dropna(subset=['Degree', 'Mutation_Count', 'Is_Metal_Binding'], inplace=True)

# Separate the data
metal_binding = data[data['Is_Metal_Binding'] == True]
non_metal_binding = data[data['Is_Metal_Binding'] == False]

# Calculate correlations and p-values
correlation_metal, p_value_metal = stats.pearsonr(metal_binding['Degree'], metal_binding['Mutation_Count'])
correlation_non_metal, p_value_non_metal = stats.pearsonr(non_metal_binding['Degree'], non_metal_binding['Mutation_Count'])

# Create output directory if it doesn't exist
output_dir = ".../"
os.makedirs(output_dir, exist_ok=True)

# Create a figure with 2 subplots, adjusted for width
fig, axs = plt.subplots(1, 2, figsize=(10, 6))  # Increased width

# Scatter plot for metal-binding proteins
sns.regplot(x='Degree', y='Mutation_Count', data=metal_binding, ax=axs[0],
            scatter_kws={'color': 'purple'}, line_kws={'color': 'purple'})
axs[0].set_xlabel('Degree', fontsize=14, fontweight='bold')
axs[0].set_ylabel('Mutation Count', fontsize=14, fontweight='bold')
axs[0].text(0.5, 0.85, f'Correlation: {correlation_metal:.2f}\nP-value: {p_value_metal:.2e}',
             fontsize=12, color='purple', fontweight='bold', ha='center', transform=axs[0].transAxes)

# Scatter plot for non-metal-binding proteins
sns.regplot(x='Degree', y='Mutation_Count', data=non_metal_binding, ax=axs[1],
            scatter_kws={'color': 'orange'}, line_kws={'color': 'orange'})
axs[1].set_xlabel('Degree', fontsize=14, fontweight='bold')
axs[1].set_ylabel('')  # Removed ylabel for the second plot
axs[1].text(0.5, 0.85, f'Correlation: {correlation_non_metal:.2f}\nP-value: {p_value_non_metal:.2e}',
             fontsize=12, color='orange', fontweight='bold', ha='center', transform=axs[1].transAxes)

# Add shared axis titles
fig.text(0.5, 0.04, 'Degree', ha='center', fontsize=14, fontweight='bold')
fig.text(0.04, 0.5, 'Mutation Count', va='center', rotation='vertical', fontsize=14, fontweight='bold')

# Add legends
axs[0].legend(['Metal Binding'], loc='upper left', fontsize=12, frameon=False)
axs[1].legend(['Non-Metal Binding'], loc='upper left', fontsize=12, frameon=False)

# Bold-face ticks and increase font size
for ax in axs:
    ax.tick_params(axis='both', which='major', labelsize=12, labelcolor='black', width=2)
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)

# Adjust layout
plt.tight_layout()
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.3)
plt.savefig(os.path.join(output_dir, 'separate_correlation_analysis_with_shared_axes_final.png'))
plt.close()

print("Final figure saved successfully!")
