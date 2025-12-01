## Five–interaction circuit screen

This folder implements the combinatorial enumeration and filtering of all
interaction circuits with exactly five terms, as described in the Methods.

### Enumeration of five–interaction topologies

- `notebooks/01_enumerate_5_interactions.ipynb` constructs all allocations of 5 interaction terms between:
  - production of \(X\) (`prod x`),
  - removal of \(X\) (`rem x`),
  - production of \(Y\) (`prod y`),
  - removal of \(Y\) (`rem y`),
  subject to the allowed monomials \(1, X, Y, X^2, Y^2, XY\).
- The total number of combinatorial possibilities is **8568**, matching:
  - the counts printed in the notebook (`number combinations` column),
  - the first entry of the filtration vector `lengths2` in cell 7.

The resulting table is stored in `data/df_5_interactions.csv`, with a `remarks`
column that documents how each combination is classified.

### “Natural” five–interaction circuits

The following filters reproduce the “natural” circuits described in the
Methods:

1. **Two‑variable dynamics**  
   Remove circuits where all 5 interactions act on a single variable
   (either X or Y), i.e. “one variable dynamics”.  
   Count after this step: **8316**.

2. **Presence of removal in both variables**  
   Remove circuits with no removal term in either X or Y.  
   Count after this step: **6990**.

3. **Presence of production in both variables**  
   Remove circuits with no production term in either X or Y.  
   Count after this step: **6210**.

These counts appear explicitly in `nb1.ipynb` (cells 1–4) and in the
filtration circles saved to `figures/venn_diagram_initial.png`.

### Additional filters toward excitable circuits

The subsequent filters reduce the “natural” set to the excitable classes:

4. **Remove X–Y symmetric duplicates**  
   In `nb1.ipynb` (cell 2), pairs of symmetric combinations are tagged
   with the remark `symmetric combination` and removed, leaving **3105**
   non‑symmetric circuits.

5. **Connectivity, feedback and phase‑space constraints**  
   The notebooks `notebooks/03_filter_by_dynamics_row1.ipynb`–`notebooks/03_filter_by_dynamics_row8.ipynb` and the merged table
   `data/df_5_interactions.csv` apply additional remarks:
   - `not connected` (disconnected circuits),
   - `no feedback` (no feedback loops),
   - `zero behaviour` / `no line of fixed points`,
   - stability‑related remarks such as `0 is not stable`,
     `another stable fp at all parameter values`.
   These filters correspond to the cumulative counts shown in `notebooks/01_enumerate_5_interactions.ipynb`
   (`lengths2`/`lengths_all`, cells 7–8) and in the venn diagrams saved
   under `figures/venn_diagram_*.png`.

6. **Stable OFF state and uniqueness of attractor**  
   For each remaining circuit, the origin \((X,Y)=(0,0)\) is required to be
   stable (or at most weakly stable) and no other stable fixed point may
   exist for all parameter values. The helper functions in
   `scripts/conditions_five.py` implement linear and center‑manifold stability
   checks that are used in the row notebooks to annotate circuits with the
   corresponding remarks.

7. **Excitability and robustness**  
   A small subset of circuits supports excitable transients away from the
   OFF state for a **robust range of parameters**. These are the circuits
   plotted in:
   - `figures/combos_remain_zoom_in_in_zero.png`,
   - `figures/circuits_remain.png`,
   and form the three excitable classes highlighted in the main figures.

Together, these scripts and notebooks implement the full pipeline from the
8568 combinatorial five–interaction circuits to the final, robustly
excitable topologies.


