# Cytokine Circuits Scan Package

This package implements the systematic scan of Cytokine-Cytokine circuits. It is structured to follow the pipeline: Generation -> Screening -> Comparison -> Feedback Variants.

## Directory Structure

- `cytokine_circuits/`: Main Python package.
    - `core/`: Core classes for Circuit definition, Integration, Plotting, and **Parameters**.
    - `generation/`, `screening/`, `comparison/`, `feedback/`: Logic modules.
- `notebooks/`: Jupyter notebooks for the analysis pipeline (**Primary Workflow**).
- `scripts/`: (Legacy/Alternative) Executable scripts.
- `data/`: Directory for input/output CSVs.

## Usage (Notebooks)

Run the notebooks in the `notebooks/` directory in order:

1.  **`01_generation.ipynb`**
    - Generates 3583 circuits (3 interactions). 
    - Output: `data/intermediate/df_all_3_interactions.csv`

2.  **`02_screening.ipynb`**
    - Screens circuits through a multi-step pipeline:
        - **Connectedness**: Ensures x and y interact.
        - **Symmetry**: Removes x-y topological duplicates.
        - **Feedback**: Ensures mutual regulation (x<->y).
        - **Divergence**: Checks stability (filters unbounded growth).
        - **Fixed Points**: Selects circuits with specific fixed point properties (single stable state).
        - **Excitability**: Checks for N-shaped nullclines (threshold behavior).
    - Classifies excitable circuits into Model Types **A** and **B**.
    - Output: `data/intermediate/df_screened_topology.csv`

3.  **`03_comparison.ipynb`**
    - Simulates screened circuits (Types A and B).
    - Computes metrics (Response Strength/Speed).
    - Performs Pareto analysis.
    - Output: `data/intermediate/df_pareto.csv`
    - *Note: Requires implementing real parameter logic in `cytokine_circuits/core/parameters.py`*

4.  **`04_feedback_analysis.ipynb`**
    - Generates **Positive (+) and Negative (-)** feedback variants for Type A and B circuits.
    - Simulates and compares them.
    - Output: `data/intermediate/df_feedback_variants.csv`

## Dependencies
- numpy, pandas, scipy, matplotlib, sympy, jupyter
