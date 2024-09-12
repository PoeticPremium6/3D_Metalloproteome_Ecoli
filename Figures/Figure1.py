###############################
#Let's start making Figure 1
#We can start with a piechart
# Load the data from the CSV files
binding_df = pd.read_csv("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv")
non_metal_df = pd.read_csv("C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\non_metal_data.csv")

# Find unique PDB_names in each DataFrame
unique_binding = binding_df['PDB_name'].unique()
unique_non_metal = non_metal_df['PDB_name'].unique()

# Count the number of unique PDB_names
num_unique_binding = len(unique_binding)
num_unique_non_metal = len(unique_non_metal)

# Print out the counts
print("Unique PDB_names in binding_data_New.csv:", num_unique_binding)
print("Unique PDB_names in non_metal_data.csv:", num_unique_non_metal)

import matplotlib.pyplot as plt

# Provided counts of unique PDB_names
num_unique_binding = 2835
num_unique_non_metal = 4506

# Data to plot
labels = 'Metal-binding', 'Metal-free'
sizes = [num_unique_binding, num_unique_non_metal]
colors = ['#7f7fff', '#ff7f7f']  # Muted shades of blue and red

# Function to format the autopct with actual value
def func(pct, allvals):
    absolute = int(pct/100.*sum(allvals))
    return "{:.1f}%\n({:d})".format(pct, absolute)

# Customize font
font = {'weight': 'bold', 'size': 30}

# Plot
plt.figure(figsize=(8, 6))  # Optional: specify the figure size
plt.pie(sizes, labels=labels, colors=colors, autopct=lambda pct: func(pct, sizes),
        startangle=0, textprops=font)
plt.axis('equal')  # Equal aspect ratio ensures that pie chart is drawn as a circle.

# Save the figure to your specified directory
figure_path = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure1\\PiechartA.png"
plt.savefig(figure_path, bbox_inches='tight')

# Show the plot
plt.show()

# Next let's show the distribution of metal-binding proteins
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Path to your CSV file
csv_file_path = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv"

# Read the CSV data into a DataFrame, set low_memory=False to prevent DtypeWarning
binding_df = pd.read_csv(csv_file_path, low_memory=False)

# Drop duplicate PDB_name and Metal_type pairs
unique_metal_types = binding_df.drop_duplicates(subset=['PDB_name', 'Metal_type'])

# Count each Metal_type
metal_type_counts = unique_metal_types['Metal_type'].value_counts()

# Total count of all metals
total_metal_count = metal_type_counts.sum()

# Select the top 10 metal types
top_10_metal_types = metal_type_counts.head(10)
other_metals_count = total_metal_count - top_10_metal_types.sum()

# Prepare data for the pie chart
sizes = list(top_10_metal_types.values) + [other_metals_count]
labels = list(top_10_metal_types.index) + ['Misc. Metals']
colors = plt.cm.tab20(np.linspace(0, 1, len(labels)))

# Create the pie chart
fig, ax = plt.subplots(figsize=(15, 8))

# Only wedges and texts are returned; no labels are added to the slices
wedges, texts = ax.pie(
    sizes,
    colors=colors,
    startangle=90
)

# Add arrows and labels for each slice outside of the pie
for i, (wedge, label) in enumerate(zip(wedges, labels)):
    angle = (wedge.theta2 - wedge.theta1) / 2.0 + wedge.theta1
    x = np.cos(np.deg2rad(angle)) * 1.2  # Adjusted for better placement
    y = np.sin(np.deg2rad(angle)) * 1.2
    ax.annotate(f'{label} ({sizes[i] / total_metal_count * 100:.1f}%)',
                xy=(x, y), xytext=(x * 1.3, y * 1.3),
                arrowprops=dict(facecolor='black', arrowstyle="->"),
                ha='center', va='center', fontsize=10, fontweight='bold')

plt.axis('equal')  # Equal aspect ratio ensures that pie chart is drawn as a circle.

# Save the figure
output_figure_path = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure1\\PiechartB.png"
plt.savefig(output_figure_path, bbox_inches='tight')
