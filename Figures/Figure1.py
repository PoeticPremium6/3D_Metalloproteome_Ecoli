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


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
file_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\binding_data_New.csv'
data = pd.read_csv(file_path)

# Count the total metal binding sites per metal
metal_counts = data['Metal_type'].value_counts()

# Setting up the color palette
# Setting up the color palette
colors = sns.color_palette('viridis', len(metal_counts))  # Replace 'viridis' with your preferred palette

# Plotting the distribution as a bar plot with a log scale
plt.figure(figsize=(12, 8))
metal_counts.plot(kind='bar', color=colors, logy=True, width=0.8)  # Adjust the width as needed
plt.title('Log-Scaled Distribution of Metal Binding Sites by Metal Type', fontsize=16, fontweight='bold')
plt.xlabel('Metal Type', fontsize=16, fontweight='bold')
plt.ylabel('Number of Binding Sites (Log Scale)', fontsize=16, fontweight='bold')
plt.xticks(rotation=45, fontsize=14, fontweight='bold')
plt.yticks(fontsize=14, fontweight='bold')
# Save the log-scaled bar plot
log_distribution_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure1\\B_metal_binding_sites_log_distribution.png'  # Adjust this path as necessary
plt.savefig(log_distribution_path)
plt.close()

print(f"Log-scaled distribution bar plot saved to: {log_distribution_path}")

#Now let's count up the EC number classification by metal type for our dataset:
#ec_class data is a basic mapping file that is readily obtained online
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Load the datasets
metals_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\refined_integrated_dataset_metals_final.csv'
non_metals_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\refined_integrated_dataset_nonmetal.csv'
ec_class_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\EC_Class.csv'

metals_data = pd.read_csv(metals_data_path)
non_metals_data = pd.read_csv(non_metals_data_path)
ec_class_data = pd.read_csv(ec_class_path)

# Label the data
metals_data['Type'] = 'Metal'
non_metals_data['Type'] = 'Non-Metal'
non_metals_data['Metal_type'] = 'Non-Metal'  # Add 'Metal_type' column for consistency

# Combine the datasets
combined_data = pd.concat([metals_data, non_metals_data], ignore_index=True)

# Extracting the first part of the EC number and format it to match EC_Class
combined_data['EC_Number_Modified'] = combined_data['EC_Number'].apply(
    lambda x: str(x).split('.')[0] + '.-.-.-' if pd.notnull(x) else None)

# Determine the top-4 metals
top_metal_types = combined_data[combined_data['Type'] == 'Metal']['Metal_type'].value_counts().head().index

# Filter the data to include only top-4 metals and non-metals
filtered_data = combined_data[
    combined_data['Metal_type'].isin(top_metal_types) | (combined_data['Type'] == 'Non-Metal')]

# Mapping modified EC numbers to their classes
mapped_data = pd.merge(filtered_data, ec_class_data, left_on='EC_Number_Modified', right_on='EC_Number', how='left')

# Data Analysis
grouped_data = mapped_data.groupby(['Metal_type', 'Class', 'Class_Function']).size().reset_index(name='Count')

# Plotting
plt.figure(figsize=(15, 10))

# Define colors
non_metal_color = 'gray'  # Color for non-metals
# Adjusting the range to start from a brighter part of the viridis colormap
metal_colors = plt.cm.viridis(np.linspace(0.2, 1, len(top_metal_types)))
metal_color_mapping = {}  # Dictionary to store metal-color mapping

# Initialize a dictionary for cumulative counts
cumulative_counts = dict.fromkeys(grouped_data['Metal_type'].unique(), 0)

# Create a stacked bar plot
for metal in grouped_data['Metal_type'].unique():
    metal_data = grouped_data[grouped_data['Metal_type'] == metal]
    color = non_metal_color if metal == 'Non-Metal' else metal_colors[top_metal_types.tolist().index(metal)]

    plt.bar(metal_data['Class_Function'] + ' (' + metal_data['Class'].astype(str) + ')',
            metal_data['Count'],
            label=metal,
            color=color,
            bottom=[cumulative_counts[m] for m in metal_data['Metal_type']])

    # Update the cumulative counts
    for m, c in zip(metal_data['Metal_type'], metal_data['Count']):
        cumulative_counts[m] += c

# Print metal-color mapping
print("Metal-Color Mapping:")
for metal, color in metal_color_mapping.items():
    print(f"{metal}: {color}")

# Enhancing Font Size and Bold Labels
plt.title('', fontsize=20, fontweight='bold')
plt.xlabel('EC Class (Function)', fontsize=20, fontweight='bold')
plt.ylabel('Counts', fontsize=20, fontweight='bold')
plt.yticks(rotation=0, fontsize=22, fontweight='bold')
plt.xticks(rotation=45, fontsize=24, fontweight='bold', ha='right')
plt.yscale('log')  # Set y-axis to logarithmic scale
plt.legend(title='Types', fontsize=16)
plt.tight_layout()

# Saving the plot
figures_directory = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure1\\'
if not os.path.exists(figures_directory):
    os.makedirs(figures_directory)

plt.savefig(os.path.join(figures_directory, 'C_EC_Classes_Top4_Metals_NonMetals.png'))

#Let's plot a description of the top-10 GO term Descriptions for each category for metals vs. non-metal
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
from matplotlib.ticker import MaxNLocator

# Load the datasets for metals and non-metals
metals_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\merged_metals_data.csv'
non_metals_data_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\merged_non_metals_data.csv'

# Read the data
metals_data = pd.read_csv(metals_data_path)
non_metals_data = pd.read_csv(non_metals_data_path)

# Function to get the top n GO terms for a given category and data
def get_top_n_go_terms(data, category, n=5):
    top_terms = (
        data[data['Category'] == category]
            .groupby('Description')
            .size()
            .reset_index(name='Count')
            .sort_values(by='Count', ascending=False)
            .head(n)
    )
    return top_terms

# Define the categories you want to plot
categories = ['molecular_function', 'biological_process', 'cellular_component']

# Function to wrap text
def wrap_labels(ax, width, axis='y'):
    if axis == 'y':
        labels = ax.get_yticklabels()
    elif axis == 'x':
        labels = ax.get_xticklabels()

    wrapped_labels = [textwrap.fill(label.get_text(), width) for label in labels]

    if axis == 'y':
        ax.set_yticklabels(wrapped_labels)
    elif axis == 'x':
        ax.set_xticklabels(wrapped_labels)

# Start plotting
fig, axes = plt.subplots(3, 2, figsize=(17, 18))  # Increased vertical size
plt.subplots_adjust(hspace=1.5, wspace=0.8)  # Adjusted vertical and horizontal space

for i, category in enumerate(categories):
    # Get the top 5 GO term descriptions for metal-binding genes
    top_metal_terms = get_top_n_go_terms(metals_data, category, n=5)

    # Get the top 5 GO term descriptions for non-metal-binding genes
    top_non_metal_terms = get_top_n_go_terms(non_metals_data, category, n=5)

    # Plot for metal-binding genes
    ax1 = sns.barplot(y='Description', x='Count', data=top_metal_terms, ax=axes[i, 0], color='blue')
    ax1.set_title(category, fontweight='bold', fontsize=15)
    ax1.set_xscale('log')
    ax1.set_xlabel('Count (Log Scale)', fontweight='bold', fontsize=12)
    ax1.set_ylabel('', fontweight='bold', fontsize=12)
    ax1.set_xticks([1, 10, 100, 1000])
    ax1.get_xaxis().set_major_formatter(plt.ScalarFormatter(useMathText=True))

    # Set bold font properties for x-tick labels
    for tick in ax1.get_xticklabels():
        tick.set_fontweight('bold')
        tick.set_fontsize(12)

    # Set bold font properties for y-tick labels
    for tick in ax1.get_yticklabels():
        tick.set_fontweight('bold')
        tick.set_fontsize(18)

    # Set tick parameters for bold tick marks
    ax1.tick_params(axis='x', which='both', width=1.5)
    ax1.tick_params(axis='y', which='both', width=1.5)

    # Wrap y-tick labels
    wrap_labels(ax1, width=26, axis='y')  # Adjust width as necessary

    # Reduce the number of y-ticks to prevent overlap
    ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))

    # Plot for non-metal-binding genes
    ax2 = sns.barplot(y='Description', x='Count', data=top_non_metal_terms, ax=axes[i, 1], color='orange')
    ax2.set_title(category, fontweight='bold', fontsize=15)
    ax2.set_xscale('log')
    ax2.set_xlabel('Count (Log Scale)', fontweight='bold', fontsize=12)
    ax2.set_ylabel('', fontweight='bold', fontsize=12)
    ax2.set_xticks([1, 10, 100, 1000])
    ax2.get_xaxis().set_major_formatter(plt.ScalarFormatter(useMathText=True))

    # Set bold font properties for x-tick labels
    for tick in ax2.get_xticklabels():
        tick.set_fontweight('bold')
        tick.set_fontsize(12)

    # Set bold font properties for y-tick labels
    for tick in ax2.get_yticklabels():
        tick.set_fontweight('bold')
        tick.set_fontsize(18)

    # Set tick parameters for bold tick marks
    ax2.tick_params(axis='x', which='both', width=1.5)
    ax2.tick_params(axis='y', which='both', width=1.5)

    # Wrap y-tick labels
    wrap_labels(ax2, width=16, axis='y')  # Adjust width as necessary

    # Reduce the number of y-ticks to prevent overlap
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=5))

# Adjust layout with padding to prevent clipping of labels
plt.tight_layout(pad=3.0, rect=[0.03, 0.03, 1, 1])  # Added padding to the left side
plt.show()

# Save the figure
output_figure_path = 'C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\Figures\\5.0\\Figure1\\E_split_top_go_terms.png'
plt.savefig(output_figure_path)
