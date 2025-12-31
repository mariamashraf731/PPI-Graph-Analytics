import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from IPython.display import display
import plotly.graph_objects as go
import random
import math
import zipfile
import os
import requests

def find_acyclic_shortest_paths(graph, source, target):
    """
    Finds the most confident (highest probability) acyclic paths between two proteins in a directed graph.

    Args:
        graph (networkx.DiGraph): The PPI network.
        source (str): UniProt ID of the source protein.
        target (str): UniProt ID of the target protein.

    Returns:
        list of tuples: Each tuple contains (path, total_path_score, interaction_weights).
    """
    try:
        # Convert probabilities to negative log values for path minimization
        for u, v, d in graph.edges(data=True):
            if d['weight'] == 0:
                d['cost'] = -math.log(d['weight']+1e-10)
            else:
                d['cost'] = -math.log(d['weight'])  # Transform confidence to cost

        # Find all shortest paths based on weight using Dijkstra's algorithm
        paths = list(nx.all_shortest_paths(graph, source=source, target=target, weight='cost', method='dijkstra'))
        #paths = list(nx.all_shortest_paths(graph, source=source, target=target, weight='weight', method='dijkstra'))
        path_details = []
        
        for path in paths:
            # Calculate the total path score (sum of edge weights)
            path_cost = math.prod(graph[path[i]][path[i+1]]['weight'] for i in range(len(path) - 1))
            # Get individual interaction weights along the path
            interaction_cost = [graph[path[i]][path[i+1]]['weight'] for i in range(len(path) - 1)]
            path_details.append((path, path_cost, interaction_cost))        
        return path_details

    except nx.NetworkXNoPath:
        # Handle the case when no path exists between the given nodes
        print(f"No path found between {source} and {target}")
        return []

# Function to write the path details to a text file
def write_paths_to_file(paths, output_file):
    """
    Writes the computed shortest paths and their details to a text file.

    Args:
        paths (list): List of paths with scores and Cost.
        output_file (str): Name of the output text file.
    """
    with open(output_file, 'w') as f:
        for path, total_score, Cost in paths:
            f.write(f"Path: {' -> '.join(path)}\n")  # Write the path in readable format
            f.write(f"Total Score: {total_score}\n")  # Write the total path score
            f.write(f"Confidence: {Cost}\n")  # Write the individual Cost of interactions
            f.write("\n")  # Add spacing between paths

# Function to visualize the sub-network formed by the shortest paths
def visualize_subnetwork(graph, paths, output_image, source, target):
    """
    Extracts and visualizes the sub-network formed by the shortest paths.

    Args:
        graph (networkx.DiGraph): The PPI network.
        paths (list): List of shortest paths.
        output_image (str): File name to save the visualization.
        source (str): Source protein ID.
        target (str): Target protein ID.
    """
    subgraph_nodes = set()  # Set to collect unique nodes from paths
    for path, _, _ in paths:
        subgraph_nodes.update(path)  # Add all nodes from each path
    
    subgraph = graph.subgraph(subgraph_nodes)  # Create a subgraph with selected nodes
    
    # Define layout for visualization
    pos = nx.spring_layout(subgraph)  
    # Create edge labels with interaction weights
    edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in subgraph.edges(data=True)}
    
    # Color nodes differently for source and target
    node_colors = ['red' if node == source else 'green' if node == target else 'lightblue' for node in subgraph.nodes]
    
    # Create the visualization
    plt.figure(figsize=(10, 8))
    nx.draw(subgraph, pos, with_labels=True, node_color=node_colors, edge_color='gray', 
            node_size=2000, font_size=10)
    nx.draw_networkx_edge_labels(subgraph, pos, edge_labels=edge_labels)
    
    plt.title("Sub-Network Visualization")
    plt.savefig(output_image)  # Save the plot as an image
    plt.show()

###################################### Data Preparation ###############################################
# Load the dataset
data_file = 'PathLinker_2018_human-ppi-weighted-cap0_75.txt'  # Update with your file path
ppi_data = pd.read_csv(data_file, sep='\t', header=1)

# Assign meaningful column names
ppi_data.columns = ['Source', 'Target', 'Edge_Weight', 'Edge_Type']
#'''
# Show Data Overview (Old Method)
print("Dataset Overview:")
display(ppi_data.head())  # Display the first few rows of the dataset
total_samples = ppi_data.shape[0]  # Total number of rows
print(f"\nTotal Number of Samples in the Dataset: {total_samples}")

# Check for duplicates
duplicates = ppi_data.duplicated()
duplicate_count = duplicates.sum()

# Display duplicate count
duplicate_data = pd.DataFrame({"Metric": ["Number of Duplicate Rows"], "Value": [duplicate_count]})
display(duplicate_data)

# If duplicates exist, display the rows
if duplicate_count > 0:
    print("\nDuplicate Rows:")
    display(ppi_data[duplicates])

# Check for missing values
missing_values = ppi_data.isnull().sum().reset_index()
missing_values.columns = ['Column', 'Missing Values']

# Display missing values as a table
print("\nMissing Values in Dataset:")
display(missing_values.style.background_gradient(cmap="Reds"))

# Check unique values in 'Source' and 'Target' columns
unique_source = set(ppi_data['Source'].unique())  # Unique tail nodes
unique_target = set(ppi_data['Target'].unique())  # Unique head nodes

# Find common nodes between Source and Target
common_nodes = unique_source.intersection(unique_target)  # Common nodes
common_nodes_count = len(common_nodes)

# Display unique counts and common nodes as a final table
final_counts = pd.DataFrame({
    "Metric": ["Unique Tail Nodes (Source)", "Unique Head Nodes (Target)", "Common Nodes"],
    "Count": [len(unique_source), len(unique_target), common_nodes_count]
})
display(final_counts)

# Display a sample of common nodes
if common_nodes_count > 0:
    common_nodes_sample = pd.DataFrame({"Common Nodes (Sample)": list(common_nodes)[:10]})
    print("\nSample Common Nodes:")
    display(common_nodes_sample)
else:
    print("\nNo common nodes between Source and Target.")

# Histogram for interaction counts
# Count the occurrences of each protein in both Source and Target columns
interaction_counts = pd.concat([ppi_data['Source'], ppi_data['Target']]).value_counts()

# Separate the top 10 and bottom 10 proteins
top_10_interactions = interaction_counts.head(10)
bottom_10_interactions = interaction_counts.tail(10)

# Plot two histograms as subplots
plt.figure(figsize=(14, 8))

# Top 10 interactions
plt.subplot(1, 2, 1)
sns.barplot(x=top_10_interactions.values, y=top_10_interactions.index, palette="Blues_d")
plt.title('Top 10 Proteins by Interaction Count', fontsize=16)
plt.xlabel('Interaction Count', fontsize=14)
plt.ylabel('Protein', fontsize=14)

# Bottom 10 interactions
plt.subplot(1, 2, 2)
sns.barplot(x=bottom_10_interactions.values, y=bottom_10_interactions.index, palette="Reds_d")
plt.title('Bottom 10 Proteins by Interaction Count', fontsize=16)
plt.xlabel('Interaction Count', fontsize=14)
plt.ylabel('Protein', fontsize=14)

# Adjust layout and show the plot
plt.tight_layout()
plt.show()


# # Find unique and uncommon proteins
# unique_source = set(ppi_data['Source'].unique())  # Unique tail nodes
# unique_target = set(ppi_data['Target'].unique())  # Unique head nodes

# Uncommon proteins
source_only = unique_source - unique_target  # Proteins exclusive to Source
target_only = unique_target - unique_source  # Proteins exclusive to Target

# Count Tail-only and Head-only Proteins
tail_only_count = len(source_only)
head_only_count = len(target_only)

# Display the counts
print(f"Number of Tail-only Proteins: {tail_only_count}")
print(f"Number of Head-only Proteins: {head_only_count}")

# Combine and label the results
uncommon_proteins = pd.DataFrame({
    "Protein": list(source_only) + list(target_only),
    "Exclusivity": ["Tail-only"] * len(source_only) + ["Head-only"] * len(target_only)
})

# Display the first few rows of uncommon proteins
print("\nUncommon Proteins:")
from IPython.display import display
display(uncommon_proteins)  # Show a sample of 20 uncommon proteins
#'''
########################################################### Graph Construction ###################################################

# Create a directed graph
G = nx.DiGraph()

# Add edges with weights
for _, row in ppi_data.iterrows():
    G.add_edge(row['Source'], row['Target'], weight=row['Edge_Weight'])

#'''
# Sample edges and check nodes
# Randomly sample 10 edges from the graph
sample_edges = random.sample(list(G.edges(data=True)), 10)
sample_nodes = set([edge[0] for edge in sample_edges] + [edge[1] for edge in sample_edges])  # Collect nodes

# Check if nodes are in the original dataset
node_checks = {node: (node in ppi_data['Source'].values or node in ppi_data['Target'].values) for node in sample_nodes}

# Prepare data for visualization
# Generate positions for nodes
pos = nx.spring_layout(G, seed=42)

# Extract edge data for visualization
edge_x, edge_y, edge_text = [], [], []
for edge in sample_edges:
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    edge_x.extend([x0, x1, None])
    edge_y.extend([y0, y1, None])
    edge_text.append(f"{edge[0]} → {edge[1]}<br>Weight: {edge[2]['weight']}")

# Create edge traces for Plotly
edge_trace = go.Scatter(
    x=edge_x,
    y=edge_y,
    line=dict(width=1.5, color='#888'),
    hoverinfo='text',
    mode='lines',
    text=edge_text
)

# Extract node data
node_x, node_y, node_text = [], [], []
for node in sample_nodes:
    x, y = pos[node]
    node_x.append(x)
    node_y.append(y)
    node_text.append(f"{node}<br>In Data: {node_checks[node]}")

# Create node traces for Plotly
node_trace = go.Scatter(
    x=node_x,
    y=node_y,
    mode='markers+text',
    text=node_text,
    hoverinfo='text',
    marker=dict(
        size=10,
        color='#1f77b4'
    )
)

# Visualize with Plotly
fig = go.Figure(
    data=[edge_trace, node_trace],
    layout=go.Layout(
        title="Sample Edges Visualization",
        titlefont_size=16,
        showlegend=False,
        hovermode='closest',
        margin=dict(b=0, l=0, r=0, t=40),
        xaxis=dict(showgrid=False, zeroline=False),
        yaxis=dict(showgrid=False, zeroline=False)
    )
)

fig.show()


# # Sample edges and check nodes
# # Randomly sample 10 edges from the graph
# sample_edges = random.sample(list(G.edges(data=True)), 10)
# sample_nodes = set([edge[0] for edge in sample_edges] + [edge[1] for edge in sample_edges])  # Collect nodes

# Check if nodes are in the original dataset
node_checks = {node: (node in ppi_data['Source'].values or node in ppi_data['Target'].values) for node in sample_nodes}

# Check for reverse edges in the original dataset
reverse_edge_checks = []
for edge in sample_edges:
    tail, head = edge[0], edge[1]
    reverse_exists = ((ppi_data['Source'] == head) & (ppi_data['Target'] == tail)).any()
    reverse_edge_checks.append((tail, head, reverse_exists))

# Prepare data for visualization
# Generate positions for nodes
pos = nx.spring_layout(G, seed=42)

# Extract edge data for visualization
edge_x, edge_y, edge_text = [], [], []
edge_weight_x, edge_weight_y, edge_weight_text = [], [], []  # To display edge weights on the plot
for edge, reverse_exists in zip(sample_edges, reverse_edge_checks):
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    edge_x.extend([x0, x1, None])
    edge_y.extend([y0, y1, None])
    edge_text.append(f"{edge[0]} → {edge[1]}<br>Weight: {edge[2]['weight']}<br>Reverse Exists: {reverse_exists[2]}")
    
    # Compute midpoints for weight labels
    edge_weight_x.append((x0 + x1) / 2)
    edge_weight_y.append((y0 + y1) / 2)
    edge_weight_text.append(f"{edge[2]['weight']:.2f}")

# Create edge traces for Plotly
edge_trace = go.Scatter(
    x=edge_x,
    y=edge_y,
    line=dict(width=1.5, color='#888'),
    hoverinfo='text',
    mode='lines',
    text=edge_text
)

# Create edge weight traces
weight_trace = go.Scatter(
    x=edge_weight_x,
    y=edge_weight_y,
    mode='text',
    text=edge_weight_text,
    hoverinfo='none',
    textfont=dict(
        size=10,
        color='red'
    )
)

# Extract node data
node_x, node_y, node_text = [], [], []
for node in sample_nodes:
    x, y = pos[node]
    node_x.append(x)
    node_y.append(y)
    node_text.append(f"{node}<br>In Data: {node_checks[node]}")

# Create node traces for Plotly
node_trace = go.Scatter(
    x=node_x,
    y=node_y,
    mode='markers+text',
    text=node_text,
    hoverinfo='text',
    marker=dict(
        size=10,
        color='#1f77b4'
    )
)

# Visualize with Plotly
fig = go.Figure(
    data=[edge_trace, weight_trace, node_trace],
    layout=go.Layout(
        title="Sample Edges Visualization with Edge Weights",
        titlefont_size=16,
        showlegend=False,
        hovermode='closest',
        margin=dict(b=0, l=0, r=0, t=40),
        xaxis=dict(showgrid=False, zeroline=False),
        yaxis=dict(showgrid=False, zeroline=False)
    )
)

fig.show()


# Retrieve the corresponding edges from the original dataset
sample_edges_df = pd.DataFrame(
    [(edge[0], edge[1], edge[2]['weight']) for edge in sample_edges],
    columns=['Source', 'Target', 'Edge_Weight']
)

# Match the sampled edges with the dataset
matched_edges = ppi_data.merge(sample_edges_df, on=['Source', 'Target', 'Edge_Weight'], how='inner')

# Display the matched edges for verification
print("\nSampled Edges from Graph:")
display(sample_edges_df)

print("\nMatched Edges in Dataset:")
from IPython.display import display
display(matched_edges)

# Retrieve the number of nodes and edges
num_nodes = G.number_of_nodes()
num_edges = G.number_of_edges()

# Display the graph details
print(f"Number of Nodes in the Graph: {num_nodes}")
print(f"Number of Edges in the Graph: {num_edges}")

# # Calculate unique nodes in Source (Tail) and Target (Head)
# unique_source = set(ppi_data['Source'].unique())  # Unique tail nodes
# unique_target = set(ppi_data['Target'].unique())  # Unique head nodes

# Find common, tail-only, and head-only nodes
common_nodes = unique_source.intersection(unique_target)  # Common nodes
tail_only_nodes = unique_source - unique_target  # Nodes in tail only
head_only_nodes = unique_target - unique_source  # Nodes in head only

# Validate the union
union_nodes = unique_source | unique_target  # Union of Source and Target nodes
calculated_union = common_nodes | tail_only_nodes | head_only_nodes  # Should match union_nodes

# Print the results
print(f"Total Nodes in Union (Direct Calculation): {len(union_nodes)}")
print(f"Total Nodes in Union (Sum of Components): {len(calculated_union)}")
print(f"Common Nodes: {len(common_nodes)}")
print(f"Nodes in Tail Only: {len(tail_only_nodes)}")
print(f"Nodes in Head Only: {len(head_only_nodes)}")

# Ensure correctness
assert union_nodes == calculated_union, "The union of nodes does not match the components!"
print("The union matches the sum of its components!")


# Generate the plot
plt.figure(figsize=(15, 15))

# Generate a layout for the nodes
pos = nx.spring_layout(G, seed=42)  # Adjust the layout as needed (e.g., circular_layout, kamada_kawai_layout)

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_size=10, node_color='blue', alpha=0.7)

# Draw edges
nx.draw_networkx_edges(G, pos, alpha=0.5, width=0.2, edge_color='gray')

# Optional: Add labels for a subset of nodes (e.g., first 50 nodes for clarity)
subset_labels = {node: node for i, node in enumerate(G.nodes()) if i < 50}
nx.draw_networkx_labels(G, pos, labels=subset_labels, font_size=8, font_color='black')

# Add title and axis
plt.title("Full Network Visualization", fontsize=20)
plt.axis('off')  # Turn off the axis
plt.show()
#'''
######################################################## Shortest Path Analysis ##########################################################

ppi_network = G

available_proteins = set(ppi_network.nodes)
print("Available proteins in network (first 10):", list(available_proteins)[:10])

# Example usage: Find shortest paths and save them to a file
source_protein = "Q9HBV2"
target_protein = "Q9Y4W2"

if source_protein not in ppi_network.nodes:
    print(f"Error: Source protein {source_protein} not found in the network.")
if target_protein not in ppi_network.nodes:
    print(f"Error: Target protein {target_protein} not found in the network.")

shortest_paths = find_acyclic_shortest_paths(ppi_network, source_protein, target_protein)

write_paths_to_file(shortest_paths, "shortest_paths.txt")

# Example usage: Visualize and save the sub-network formed by shortest paths
print("Path_1")
print(shortest_paths[0])
visualize_subnetwork(ppi_network, shortest_paths, "subnetwork_enhanced.png", source_protein, target_protein)

############################################################### Connectivity Analysis & protein UniProt_ID-gene_name Conversion Map #########################################################

# Function to compute the degree of a single protein
def compute_single_protein_degree(graph, protein_id):
    return graph.degree(protein_id)

# Function to compute degrees for a set of proteins
def compute_set_of_protein_degrees(graph, protein_ids):
    return {protein_id: graph.degree(protein_id) for protein_id in protein_ids}

# Function to create a histogram of protein degrees
def create_degree_histogram(graph, output_folder):
    degrees = [deg for _, deg in graph.degree()]
    plt.hist(degrees, bins=range(1, max(degrees) + 2), edgecolor='black')
    plt.title("Protein Degree Distribution")
    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    output_path = os.path.join(output_folder, "degree_histogram.png")
    plt.savefig(output_path)
    print(f"Degree histogram saved to {output_path}")
    plt.show()

# Function to create a bar chart for the top 10 protein degrees
def create_degree_histogram_with_names(graph, output_folder):
    degrees = dict(graph.degree())
    sorted_proteins = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
    top_10_proteins = sorted_proteins[:10]
    protein_names = [protein for protein, _ in top_10_proteins]
    top_10_degrees = [degree for _, degree in top_10_proteins]
    plt.figure(figsize=(10, 6))
    plt.bar(protein_names, top_10_degrees, color='skyblue', edgecolor='black')
    plt.title("Top 10 Protein Degree Distribution")
    plt.xlabel("Protein Names")
    plt.ylabel("Degree")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    output_path = os.path.join(output_folder, "top_10_proteins.png")
    plt.savefig(output_path)
    print(f"Top 10 proteins histogram saved to {output_path}")
    plt.show()

# Function to rank proteins by degree
def rank_proteins_by_degree(graph, output_folder):
    degrees = dict(graph.degree())
    ranked = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
    output_file = os.path.join(output_folder, "ranked_degrees.txt")
    with open(output_file, 'w') as file:
        for protein, degree in ranked:
            file.write(f"{protein}\t{degree}\n")
    print(f"Ranked degrees saved to {output_file}")

# Function to list connected proteins
def list_directly_connected_proteins(graph, protein_id, output_folder):
    """
    List all directly connected proteins to a randomly selected protein along with their interaction weights.

    Parameters:
        graph (networkx.Graph): The protein interaction graph.
        protein_id (str): The randomly selected protein ID.
        output_folder (str): Directory to save the output file.
    """
    try:
        if protein_id not in graph:
            raise ValueError(f"Protein {protein_id} not found in the graph.")

        # Degree of the selected protein
        degree = graph.degree(protein_id)

        # Get all directly connected proteins
        connected_proteins = graph[protein_id]

        # Create the output file
        output_file = os.path.join(output_folder, f"{protein_id}_connected_proteins.txt")

        with open(output_file, "w") as file:
            file.write(f"Protein: {protein_id}\n")
            file.write(f"Degree: {degree}\n\n")
            file.write("Connected Proteins:\n")

            # List directly connected proteins and their weights
            for neighbor, attributes in connected_proteins.items():
                interaction_weight = attributes.get("weight", "N/A")
                file.write(f"{neighbor}\tWeight: {interaction_weight}\n")

        print(f"Directly connected proteins for {protein_id} saved to {output_file}")

    except Exception as e:
        print(f"Error: {e}")


# Function to fetch gene names from UniProt API
def fetch_gene_name(uniprot_id):
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            return data.get("genes", [{}])[0].get("geneName", {}).get("value", "Unknown")
        return "Unknown"
    except Exception as e:
        print(f"Error fetching gene name for {uniprot_id}: {e}")
        return "Unknown"

# Function to map UniProt IDs to gene names
def map_uniprot_to_gene_name(protein_ids, output_folder):
    output_file = os.path.join(output_folder, "uniprot_to_gene_map.txt")
    try:
        with open(output_file, "w") as file:
            file.write("UniProt_ID\tGene_Name\n")
            for pid in protein_ids:
                gene_name = fetch_gene_name(pid)
                file.write(f"{pid}\t{gene_name}\n")
        print(f"UniProt to Gene Name map saved to {output_file}")
    except Exception as e:
        print(f"Error saving UniProt mapping: {e}")

# Load the PPI network
extract_to_dir = "."
ppi_graph = G
print(f"Graph loaded successfully!")
print(f"Number of nodes: {ppi_graph.number_of_nodes()}")
print(f"Number of edges: {ppi_graph.number_of_edges()}")

# Top 10 proteins by degree
degrees = dict(ppi_graph.degree())
top_10_proteins = [protein for protein, _ in sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:10]]

# Randomly select a protein from the top 10
random_protein = random.choice(top_10_proteins)
print(f"Randomly selected protein: {random_protein}")

# List directly connected proteins
list_directly_connected_proteins(ppi_graph, random_protein, extract_to_dir)

# Additional analysis (if needed)
create_degree_histogram(ppi_graph, extract_to_dir)
create_degree_histogram_with_names(ppi_graph, extract_to_dir)
rank_proteins_by_degree(ppi_graph, extract_to_dir)

############################################################# Adjacency Matrix ########################################################

# Create the weighted adjacency matrix
weighted_adj_matrix = nx.to_pandas_adjacency(G, weight='weight')

# Convert to an unweighted graph and create the unweighted adjacency matrix
G_unweighted = nx.DiGraph()
for u, v in G.edges():
    G_unweighted.add_edge(u, v)  # Add edges without weights
unweighted_adj_matrix = nx.to_pandas_adjacency(G_unweighted, dtype=int)

# Save the unweighted adjacency matrix to a file
unweighted_adj_matrix.to_csv("unweighted_adjacency_matrix.csv")
print("Unweighted adjacency matrix saved to 'unweighted_adjacency_matrix.csv'.")

# Visualize the weighted adjacency matrix
plt.figure(figsize=(12, 10))
sns.heatmap(
    weighted_adj_matrix.iloc[:50, :50],  # Show the first 50 nodes for clarity
    cmap='viridis',
    cbar=True,
    square=True
)
plt.title("Weighted Adjacency Matrix (First 50 Nodes)", fontsize=16)
plt.xlabel("Target Nodes")
plt.ylabel("Source Nodes")
plt.show()

# Visualize the unweighted adjacency matrix
plt.figure(figsize=(12, 10))
sns.heatmap(
    unweighted_adj_matrix.iloc[:50, :50],  # Show the first 50 nodes for clarity
    cmap=['white', 'black'],  # Two-color map: white (0), black (1)
    cbar=False,  # Disable the color bar for binary data
    square=True
)
plt.title("Unweighted Adjacency Matrix (First 50 Nodes)", fontsize=16)
plt.xlabel("Target Nodes")
plt.ylabel("Source Nodes")
plt.show()
