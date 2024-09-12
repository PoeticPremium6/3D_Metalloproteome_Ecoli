#Let's move to Figure 2 which focuses on metals, mutations, residues, and proximity
#The plot displays the distances between metal ions and nearby residues in PDB structures
#X-Axis: Distance between metals and residues in Ångstroms.
#Y-Axis: Frequency of occurrence for each distance.
#Histogram Bars: Number of occurrences in each distance range.
#KDE Line: Smoothed representation of distance distribution.
#Key Insights:
#Peak Distances: Common interaction distances.
#Spread: Variation in interaction distances.
#Tail: Indicates residues farther from metals but still nearby.
#Average Distance: Central tendency of interactions.
#Outliers: Unusually large distances warranting investigation.
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the binding data from your local system
binding_df = pd.read_csv("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv")  # Update the path as needed

# Filter the dataframe for distances less than or equal to 5 Angstroms
binding_df_filtered = binding_df[binding_df['Distance'] <= 5]

# Count the frequency of interactions before and after 3 Å
total_interactions_before_3A = binding_df_filtered[binding_df_filtered['Distance'] <= 3].shape[0]
total_interactions_after_3A = binding_df_filtered[binding_df_filtered['Distance'] > 3].shape[0]

# Print the total frequency of interactions before and after 3 Å
print("Total interactions before 3 Å:", total_interactions_before_3A)
print("Total interactions after 3 Å:", total_interactions_after_3A)

# Set style
sns.set_style("whitegrid")

# Create the figure
plt.figure(figsize=(14, 8))

# Create the histogram with adjusted bins and a kernel density estimate
# This time, considering only data up to 5 Angstroms, spread equally
sns.histplot(binding_df_filtered['Distance'], bins=50, kde=True, color='skyblue', edgecolor='black')

# Calculate mean and median distances and add them as vertical lines for the filtered dataset
mean_distance = binding_df_filtered['Distance'].mean()
median_distance = binding_df_filtered['Distance'].median()
plt.axvline(mean_distance, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_distance:.2f} Å')
plt.axvline(median_distance, color='green', linestyle='-', linewidth=2, label=f'Median: {median_distance:.2f} Å')

# Add threshold lines for significant interaction distances within 5 Å
plt.axvline(2, color='purple', linestyle='--', linewidth=2, label='Direct Coordination Threshold (2 Å)')
plt.axvline(3, color='purple', linestyle='--', linewidth=2, label='Possible Coordination Threshold (3 Å)')
plt.axvline(5, color='orange', linestyle='--', linewidth=2, label='Interaction Cutoff (5 Å)')

# Add title and labels
#plt.title('Distribution of Distances (≤5 Å) Between Metals and Residues')
plt.ylabel('Frequency', fontsize=16, fontweight='bold')
plt.xlabel('Distance (Å)', fontsize=16, fontweight='bold')

# Adjusting the x-axis limit to make sure it's equally spread up to 5 Å
plt.xlim(0, 5)

# Customize x and y ticks
plt.xticks(fontsize=14, fontweight='bold')
plt.yticks(fontsize=14, fontweight='bold')

# Add a legend to the plot
legend = plt.legend(fontsize=14, title_fontsize=16)  # Add the legend
legend.get_title().set_fontweight('bold')  # Set the font weight for the legend title
# Set bold font properties for the legend text
for text in legend.get_texts():
    text.set_fontweight('bold')
    text.set_fontsize(14)
# Tighten the layout
plt.tight_layout()

# Save the figure
plt.savefig("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure2\\A_distance_distribution_equal_spread.png")  # Update the path as needed

# Show the plot
plt.show()
