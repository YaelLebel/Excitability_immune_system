# Release Preview: Immune System Excitability Analysis Codebase

We are excited to announce the release of the code and data repository accompanying our paper on excitable circuits in the immune system. This release offers a fully reproducible and structured implementation of our methods, from circuit enumeration to clinical data analysis.

## 🚀 Key Highlights

*   **Systematic Circuit Evaluation**: Access the complete pipeline for enumerating and screening 3,583 cytokine circuits and 8,568 cell-cell circuits to identify excitable topologies.
*   **Network Motif Discovery**: Tools to explore the cytokine and cell-cell interaction networks derived from ImmunoGlobe, identifying statistically enriched excitable motifs.
*   **COVID-19 Dynamics Modeling**: Replicate our ODE modeling of longitudinal patient data, capturing the interplay between pro-inflammatory markers and IL-10.
*   **Methodological Transparency**: The codebase is architected to mirror the paper's structure, ensuring every figure and result can be traced back to its source code.

## 📂 Repository Structure
The repository layout has been optimized for navigation, with top-level directories directly mapping to the paper's method sections:

*   `Cytokine_circuits_scan/`: Systematic enumerations (Section 1)
*   `Network_scan/`: Motif search in cytokine & cell networks (Sections 2 & 7)
*   `COVID_Cytokine_Dynamics/`: Clinical data analysis & modeling (Sections 3 & 4)
*   `Comparison_to_classical_models/`: Benchmarking (Section 5)
*   `Cell_circuits_scan/`: Cell-cell circuit analysis (Sections 6 & 8)

## 🛠 Getting Started
Please refer to the `README.md` for a detailed mapping of method sections to code directories and instructions on reproducing specific results.
`