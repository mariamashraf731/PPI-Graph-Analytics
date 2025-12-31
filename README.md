# 🕸️ PPI Graph Analytics: Protein-Protein Interaction Network Analysis

![Language](https://img.shields.io/badge/Python-3.8%2B-blue)
![Library](https://img.shields.io/badge/Graph-NetworkX-orange)
![Domain](https://img.shields.io/badge/Bioinformatics-PPI-green)
![Analysis](https://img.shields.io/badge/Method-Graph%20Theory-purple)

## 📌 Project Overview
This project applies **Graph Theory** and **Network Science** to analyze Human Protein-Protein Interaction (PPI) networks. Using **NetworkX**, it explores the topological properties of the proteome to identify essential proteins ("hubs") and discover potential signaling pathways.

The analysis provides biological insights by modeling proteins as **Nodes** and their interactions as **Edges**, allowing for the application of advanced graph algorithms.

## ⚙️ Key Analyses & Algorithms
1.  **Network Topology Analysis:**
    * Constructed a weighted undirected graph from PPI data.
    * Calculated **Degree Distribution** to identify scale-free properties.
    * Computed **Graph Density** and **Average Clustering Coefficient**.
2.  **Hub Identification (Centrality Measures):**
    * Used **Degree Centrality** to find the most connected proteins (e.g., UBC, APP).
    * Ranked proteins to highlight potential drug targets.
3.  **Pathfinding & Signaling:**
    * Implemented **Dijkstra’s Algorithm** (via NetworkX) to find **Shortest Paths** between specific proteins.
    * Analyzed connectivity for specific targets (e.g., **P04629** - NTRK1).
4.  **Visualization:**
    * Generated subnetwork graphs to visualize local interactions.
    * Plotted Degree Histograms on Log-Log scales.

## 📂 Dataset
* **Source:** PathLinker 2018 Human PPI (Weighted).
* **Mapping:** UniProt IDs mapped to Gene Symbols for readability.

## 🚀 How to Run
1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/mariamashraf731/PPI-Graph-Analytics.git](https://github.com/mariamashraf731/PPI-Graph-Analytics.git)
    ```
2.  **Install Dependencies:**
    ```bash
    pip install -r requirements.txt
    ```
3.  **Run the Analysis:**
    ```bash
    jupyter notebook notebooks/PPI_Analysis_Main.ipynb
    ```

## 📄 Documentation
For a detailed explanation of the graph metrics and biological interpretations, refer to the [Final Report](docs/PPI_Analysis_Report.pdf).
