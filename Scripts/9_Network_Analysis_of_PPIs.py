import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

base_directory = '.../'  # Adjust this path if necessary

# Load data
def load_data(file_path):
    return pd.read_csv(file_path, low_memory=False)

# Load mutation and PPI data
mutation_data = load_data(base_directory + 'categorized_dataset.csv')
ppi_data = load_data(base_directory + 'categorized_dataset.csv')  # Assuming this file contains PPI info

# Create a PPI network and a dictionary for mutation counts
G = nx.Graph()
mutation_counts = {}

# Add nodes and edges from PPI data
for _, row in ppi_data.iterrows():
    protein = row['Uniprot_ID']
    if pd.notna(protein) and pd.notna(row['Known_PPI']):
        interactions = str(row['Known_PPI']).split(';')  # Adjust delimiter if needed
        for interact in interactions:
            if interact:
                G.add_edge(protein, interact)

# Count mutations for each protein
for _, row in mutation_data.iterrows():
    protein = row['Uniprot_ID']
    mutation_counts[protein] = mutation_counts.get(protein, 0) + 1  # Increment count

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

# Normalize mutation counts and apply a logarithmic scale
normalized_counts = np.log([mutation_counts.get(node, 0) + 1 for node in G.nodes()])
sizes = (normalized_counts - min(normalized_counts) + 1) * 200  # Adjust the factor to scale sizes

# Draw the network with color coding and node sizing
plt.figure(figsize=(20, 20))
pos = nx.spring_layout(G, k=0.016*np.sqrt(len(G.nodes())), iterations=50)

# Draw nodes with inverted viridis colormap for high values dark and low values light
nodes = nx.draw_networkx_nodes(G, pos, node_color=normalized_counts, node_size=sizes, cmap='viridis_r', alpha=0.9)

# Draw edges
edges = nx.draw_networkx_edges(G, pos, alpha=0.4, width=0.8)
