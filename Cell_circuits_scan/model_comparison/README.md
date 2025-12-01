## Model comparison and excitable archetypes

This folder contains notebooks that instantiate specific excitable
archetypes and compare their dynamical properties within the general
two–population framework defined in `Base/`.

### `model_comparison.ipynb`

This notebook defines three archetype models (A1–A3) of two interacting
populations \(X\) and \(Y\) of the form
\[
\frac{dX}{dt} = X \, P_X(X,Y)\left(1-\frac{X}{C_X}\right) - X \, R_X(X,Y),
\qquad
\frac{dY}{dt} = Y \, P_Y(X,Y)\left(1-\frac{Y}{C_Y}\right) - Y \, R_Y(X,Y),
\]
with production and removal polynomials constructed from monomials up to
total degree 2, as in the Methods.

The notebook:
- builds each archetype using `DynamicModel2D` with explicit expressions
  for \(dX/dt\) and \(dY/dt\),
- computes and plots nullclines and typical trajectories
  (phase portraits),
- integrates responses to different input amplitudes and records:
  - maximum \(X\) output,
  - time to maximum \(X\),
  - maximum \(Y\) output,
  - refractory period (time until \(Y\) decays below a threshold),
  - early exponential growth rates.

The summary DataFrames
`df_integrate_properties_a*_b*.csv`, `df_comparison_a*_b*.csv` and
`df_comparison_maximize_a*_b*.csv` store these metrics and are used to
generate the sensitivity, refractory period and Pareto‑front figures in
the paper.

### `archetype_1/archetype_1.ipynb`

This notebook applies the same set of metrics to archetype 1 and to
representative excitable circuits selected from the five‑interaction
screen:

- the file imports five‑interaction circuits from
  `../five_interactions/row*.csv`,
- reconstructs the corresponding `Combo_CC` models,
- uses `Combo_ODE` to integrate trajectories across a range of input
  magnitudes,
- compares the input–output curves, sensitivity and refractory period of
  archetype 1 to those of the excitable classes (`cms_arc1`, `cms_arc2`),
  saving the resulting figures under `figures/`.

These analyses directly implement the “Screening for excitable circuit”
and “Comparison between excitable archetypes” steps described in the
Methods.


