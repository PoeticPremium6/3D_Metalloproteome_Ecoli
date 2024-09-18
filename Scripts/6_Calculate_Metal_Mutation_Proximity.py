#Let's begin to study more about metal residue proximity, but integrate our mutation data
#Results will be used for Figures 2D/E and eventually Figure 3
import csv
import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import math

# Directories for input and output data
binding_data_file = ".../binding_data_New.csv"
mutation_files = [
    ".../ALE_Mutations.csv",
    ".../LTEE_Mutations.csv"
]
output_directory = ".../"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Load metal-binding data
def load_csv(file_path, delimiter=','):
    with open(file_path, 'r') as f:
        return [row for row in csv.DictReader(f, delimiter=delimiter)]

binding_data = load_csv(binding_data_file)

# Load mutation datasets (ALE and LTEE)
mutation_data = []
for mutation_file in mutation_files:
    mutation_data.extend(load_csv(mutation_file, delimiter=';'))

# Find matching mutations
matching_mutations = [
    {"mutation": mutation, "binding": binding}
    for mutation in mutation_data
    for binding in binding_data
    if mutation['AA_residue'] == binding['Residue_number']
]

# Count mutations per metal type
mutations_per_metal = defaultdict(int)
for entry in matching_mutations:
    mutations_per_metal[entry['binding']['Metal_type']] += 1

# Save the mutations per metal data to a pickle file
pickle_file = os.path.join(output_directory, "mutations_per_metal.pkl")
with open(pickle_file, 'wb') as f:
    pickle.dump(mutations_per_metal, f)

# Function to plot the data
def plot_data(df, save_path, title, xlabel, ylabel):
    plt.figure(figsize=(12, 8))
    plt.bar(df['Metal'], df['Count'], log=True)
    plt.xlabel(xlabel, fontsize=12, fontweight='bold')
    plt.ylabel(ylabel, fontsize=12, fontweight='bold')
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xticks(rotation=90, fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')
    plt.savefig(save_path)
    plt.show()

# Read and plot mutation per metal data
df = pd.DataFrame(list(mutations_per_metal.items()), columns=['Metal', 'Count']).sort_values(by='Count', ascending=False)
plot_save_path = os.path.join(output_directory, "Count_Metal_Mutations.png")
plot_data(df, plot_save_path, 'Descending Order Metal Counts on a Log Scale', 'Metal', 'Count (log scale)')

# Function to compute distance between coordinates
def compute_distance(coord):
    x, y, z = map(float, coord.replace("(", "").replace(")", "").split(","))
    return math.sqrt(x**2 + y**2 + z**2)

# Compute distances between mutations and metal-binding sites
distances = []
for mutation in mutation_data:
    closest_distance = float('inf')
    for binding in binding_data:
        if mutation['AA_residue'] == binding['Residue_number']:
            distance = compute_distance(binding['Residue_coord'])
            if distance < closest_distance:
                closest_distance = distance
    if closest_distance != float('inf'):
        distances.append(closest_distance)

# Count mutations near metal-binding sites (distance threshold of 5Å)
def mutations_near_metal(metalloproteins, mutations, threshold=5):
    count = 0
    for metal in metalloproteins:
        metal_residue = int(metal['Residue_number'])
        for mutation in mutations:
            if abs(metal_residue - int(mutation['AA_residue'])) <= threshold:
                count += 1
    return count

ale_count = mutations_near_metal(binding_data, [m for m in mutation_data if 'ALE' in m['Dataset']])
ltee_count = mutations_near_metal(binding_data, [m for m in mutation_data if 'LTEE' in m['Dataset']])

# Save mutation counts near metal-binding sites
mutation_count_pickle = os.path.join(output_directory, "Mutations_Near_Metal_Counts.pkl")
with open(mutation_count_pickle, 'wb') as f:
    pickle.dump({'ALE': ale_count, 'LTEE': ltee_count}, f)
