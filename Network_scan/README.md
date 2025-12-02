# Network Scan Workflow

This directory now mirrors the workflow described in Methods §7 (“Network Search for Excitable Cell-Cell Motifs”). It captures every stage from parsing ImmunoGlobe interaction strings to motif enrichment testing on randomized graphs.

## Directory Overview
```
Network_scan/
  README.md
  Cytokine_network/
    data_raw/
    data_processed/
    notebooks/
  data_raw/
    immunoglobe_edges.csv          # ImmunoGlobe edge annotations (cells, cytokines, actions, processes)
    immunoglobe_nodes.csv          # Node metadata used to restrict to immune/leukocyte populations
  data_processed/
    filtered_cell_network.csv                      # Signed cell-cell edges with mediators retained
    filtered_cell_network_filtered_actions.csv     # Activate/Survive-only subset for graph flattening
  notebooks/
    01_construct_cell_cell_network.ipynb           # Methods §7.1 implementation
    02_flatten_motif_enrichment.ipynb              # Methods §7.2–§7.4 implementation
```

## Mapping to Methods §7
- **§7.1 – Construction of the Cell-Cell Interaction Network**  
  `notebooks/01_construct_cell_cell_network.ipynb` loads the raw inputs from `data_raw/*.csv`, removes ambiguous actions (e.g., Polarize), filters to immune cells, excludes non-leukocyte/tumor entries, and classifies interactions as direct or cytokine-mediated. The resulting mediator-aware, signed edges are exported to `data_processed/filtered_cell_network.csv` plus the Activate/Survive subset `data_processed/filtered_cell_network_filtered_actions.csv`.
- **§7.2 – Graph Flattening and Annotating**  
  `notebooks/02_flatten_motif_enrichment.ipynb` defines `create_flattened_graph_from_data_weights`, aggregates every Source→Target pair, sets `weight_positive` when any Activate/Survive edge exists, sets `weight_negative` when any Kill/Inhibit edge exists, and stores the underlying mediator tables for downstream inspection.
- **§7.3 – Motif Identification (LAD search)**  
  The same notebook implements `count_excitable_motif`, `find_excitable_motif`, and `find_excitable_motif_b`, which encode the LAD template (self-loop, positive forward edge, negative feedback edge). Matches are restricted to positive weights on X→X and X→Y and a negative weight on Y→X, yielding 21 excitable motifs with full mediator context.
- **§7.4 – Statistical Enrichment**  
  `rewire_graph_unique_edges` generates randomized graphs by swapping edge targets while preserving positive/negative labels and degree distributions. The notebook runs 1,000 such randomizations, recounts the LAD motif in each, and compares the observed count (21) to the null distribution (p = 0.04).

## File Inventory (Cell-Cell Workflow)
| Path | Purpose | Methods reference |
| --- | --- | --- |
| `data_raw/immunoglobe_edges.csv` | Raw ImmunoGlobe interaction records (cell/cytokine/action/process strings). | §7.1 |
| `data_raw/immunoglobe_nodes.csv` | Node metadata table to identify immune/leukocyte cell types. | §7.1 |
| `data_processed/filtered_cell_network.csv` | Filtered cell-cell edges (direct + mediated) retaining mediators and signed actions. | §7.1–§7.2 |
| `data_processed/filtered_cell_network_filtered_actions.csv` | Flattened Activate/Survive subset for building binary positive weights. | §7.2 |
| `notebooks/01_construct_cell_cell_network.ipynb` | Cleans raw data, excludes ubiquitous mediators, and exports processed CSVs. | §7.1 |
| `notebooks/02_flatten_motif_enrichment.ipynb` | Builds weighted igraphs, performs LAD motif discovery, and runs 1,000 rewired-network tests. | §7.2–§7.4 |

## Cytokine Network Subdirectory
`Cytokine_network/` mirrors the same structure for cytokine-only interactions:

```
Cytokine_network/
  data_raw/
    immunoglobe_edges.csv          # (Optional copy or symbolic link to root-level raw edges)
    immunoglobe_nodes.csv
  data_processed/
    filtered_cytokine_network.csv              # Cytokine→cytokine edges inferred via one-step mediation
  notebooks/
    01_construct_cytokine_interaction_network.ipynb
    02_flatten_motif_enrichment.ipynb
```

- `notebooks/01_construct_cytokine_interaction_network.ipynb` (previously `nwm_cyt_init.ipynb`) filters ImmunoGlobe to cytokine actors, assembles cytokine→cytokine interactions, and exports `data_processed/filtered_cytokine_network.csv`.
- `notebooks/02_flatten_motif_enrichment.ipynb` adapts the flattening + LAD motif pipeline to the cytokine network, enabling motif discovery/enrichment analyses that parallel Methods §7 but restricted to cytokines.

Use the same reproducibility checklist—run Notebook 01 to refresh processed cytokine edges, then execute Notebook 02 for motif counts and randomization tests specific to cytokine circuits.

## Reproducibility Checklist
1. **Notebook 01 (§7.1)** – Execute all cells to regenerate `filtered_cell_network*.csv`. Confirm that the action whitelist {Activate, Survive, Secrete, Inhibit, Kill} and mediator filtering (e.g., removal of ubiquitous cytokines such as IFN-γ, TGF-β) match the manuscript description.
2. **Notebook 02 (§7.2–§7.3)** – Load `data_processed/filtered_cell_network.csv`, rebuild the weighted graph via `create_flattened_graph_from_data_weights`, and verify that `find_excitable_motif` returns 21 motifs along with mediator annotations.
3. **Notebook 02 (§7.4)** – Run the `rewire_graph_unique_edges` randomization loop for 1,000 iterations, recompute motif counts for each randomized graph, and confirm the enrichment statistic (observed 21 motifs; p ≈ 0.04 versus the null distribution).

Maintaining this structure keeps the computation synchronized with Methods §7 and ensures every file’s role in the pipeline is explicit.
