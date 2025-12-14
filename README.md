# Excitability Immune System

This directory contains the code and data corresponding to the analysis of excitable circuits in the immune system. The structure of the codebase is organized to align with the methods sections of the associated paper.

## Directory to Methods Mapping

The following table outlines how each directory corresponds to the specific sections of the paper's methods:

| Paper Section | Method Description | Corresponding Directory |
| :--- | :--- | :--- |
| **Section 1** | Systematic Scan of Cytokine-Cytokine Circuits | `Cytokine_circuits_scan/` |
| **Section 2** | Network Search for Excitable Cytokine Motifs | `Network_scan/Cytokine_network/` |
| **Section 3** | Analysis of Longitudinal type I interferon and SARS-CoV-2 Patient Data | `COVID_Cytokine_Dynamics/` |
| **Section 4** | Modeling Pro-Inflammatory Markers and IL-10 Dynamics in COVID-19 Patients | `COVID_Cytokine_Dynamics/` |
| **Section 5** | Comparison of immune response length to viral infection duration | `Comparison_to_classical_models/` |
| **Section 6** | Systematic Scan of Cell-Cell Circuits | `Cell_circuits_scan/` |
| **Section 7** | Network Search for Excitable Cell-Cell Motifs | `Network_scan/Cell_network/` |
| **Section 8** | Circuit to target and sensitivity analysis | `Cell_circuits_scan/` |

Each directory contains the relevant scripts, notebooks, and data required to reproduce the analysis described in its corresponding section.

## Dependencies

Some parts of this codebase rely on the `DynamicModel_Package`. For convenience, a copy of this package is included in the `libs/` directory.

To use the code that depends on this package, ensure it is in your Python path. You can do this by setting the `PYTHONPATH` environment variable:

```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)/libs
```

Or by installing it in editable mode:

```bash
pip install -e libs/DynamicModel_Package
```

