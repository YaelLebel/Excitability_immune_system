## Overview

This folder implements the computational framework used to screen two–cell-type interaction circuits for excitable dynamics, as described in the Methods.

We consider two interacting populations \(X\) and \(Y\) with dynamics
\[
\frac{dX}{dt} = X \, P_X(X,Y)\left(1-\frac{X}{C_X}\right) - X \, R_X(X,Y),
\qquad
\frac{dY}{dt} = Y \, P_Y(X,Y)\left(1-\frac{Y}{C_Y}\right) - Y \, R_Y(X,Y),
\]
where:
- \(P_X, P_Y\) are production polynomials,
- \(R_X, R_Y\) are removal polynomials,
- \(C_X, C_Y\) are carrying capacities.

The polynomial terms implement:
- constant background proliferation/removal,
- linear and quadratic self-activation or limitation,
- cross-activation or inhibition,
- cooperative inhibition (e.g. terms in \(XY, X^2, Y^2\)).
Only monomials up to total degree 2 in \(X,Y\) are used inside the brackets, which correspond to the coefficient notation \(a_{ij} X^i Y^j\) with \(0 \le i+j \le 2\) and non‑negative coefficients.

### Folder structure

- `Base/`  
  Core classes and utilities that encode the generic two–population ODE model and the combinatorial space of interaction terms:
  - `Combo_CC_Base.py`: implements `Combo_CC`, which builds \(P_X, P_Y, R_X, R_Y\) polynomials and the full ODEs with carrying capacities and optional positive parameters (`dp0`, `dp1`, …). This enforces \(X=0\) and \(Y=0\) as absorbing states and bounds the dynamics.
  - `Combo_Base.py`, `coeff_base.py`, `Combo_utils.py`: encode the coefficient representation (`+/- a_ij`, `+/- b_ij`), generate all allowed combinations of interaction monomials and reconstruct the corresponding polynomials \(P(X,Y)\) and \(Q(X,Y)\).
  - `Combo_ode_base.py`, `Combo_plot_base.py`, `utils.py`: utilities for numerical integration, phase‑portrait plotting, drawing filtration/venn diagrams and schematic circuit diagrams.

- `five_interactions/`  
  Scripts, notebooks and data to enumerate and filter all circuits with exactly five interaction terms:
  - `notebooks/01_enumerate_5_interactions.ipynb` and `notebooks/03_filter_by_dynamics_row*.ipynb` construct all combinatorial allocations of 5 terms between production/removal of \(X\) and \(Y\) (total 8568 circuits), then successively filter:
    1. Remove “one‑variable” circuits (all 5 terms in one population) → 8316.
    2. Require at least one removal term in each population → 6990.
    3. Require at least one production term in each population → 6210.
    4. Remove circuits related by \(X \leftrightarrow Y\) symmetry → 3105.
    5. Remove disconnected circuits, circuits without feedback, circuits with unstable or non‑hyperbolic origin, and circuits where another stable fixed point exists for all parameter values.  
       These filters correspond to the additional counts shown in the notebooks and figures, and reduce the set to the final excitable classes.
  - `data/df_5_interactions.csv`, `data/row*.csv`: tabulate the combinatorial classes and the sequence of remarks/filters applied to each.
  - `scripts/conditions_five.py`: helper routines for linear and center‑manifold stability analysis of the origin, used when checking that the “OFF” state at \((0,0)\) is stable.

- `model_comparison/`  
  Notebooks that instantiate specific excitable archetypes and quantify their dynamical properties:
  - `model_comparison.ipynb`: defines three archetype models (`A1–A3`) within the same general framework and computes input–output curves, sensitivity, refractory period, and time‑to‑peak.
  - `archetype_1/archetype_1.ipynb`: applies the same metrics to archetype 1 and representative excitable circuits obtained from the five‑interaction screening. The CSVs in this folder (`df_traj_*`, `df_integrate_properties_*`, `df_comparison_*`) store the trajectories and summary statistics used in the figures.

Together, these components implement the full pipeline described in the Methods: from defining the polynomial interaction space, through enumerating and filtering all five‑interaction circuits, to detailed dynamical comparison of the excitable archetypes.


