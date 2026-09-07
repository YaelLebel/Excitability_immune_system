# Circuit to Target and Sensitivity Analysis

This directory implements the circuit-to-target analysis of the excitable
cell-cell circuit, corresponding to the Methods section **"Circuit to target and
sensitivity analysis"** and to **Fig. 6** of the paper. It asks which parameters
of the excitable circuit are effective targets for reducing the immune response
(as in autoimmunity) or enhancing it (as in cancer or hypoimmunity), and which
dynamical regimes the circuit can occupy when parameters are changed by large
amounts.

## The model

The analysis uses a reduced version of the excitable cell-cell architecture
found by the five-interaction screen in `Cell_circuits_scan/`:

```
dX/dt = X [ a X (1 - X/c) - b Y - d ]
dY/dt = Y [ f X - g ]
```

| Parameter | Meaning | Units |
| :--- | :--- | :--- |
| `a` | autocrine (self-activation) strength of the effector `X` | 1 / (X t) |
| `b` | inhibition of `X` by the inhibitor `Y` | 1 / (Y t) |
| `c` | carrying capacity of `X` | X |
| `d` | removal rate of `X` | 1 / t |
| `f` | induction of the inhibitor `Y` by `X` | 1 / (X t) |
| `g` | turnover (removal) rate of `Y` | 1 / t |

As in the general cell-cell class (`Cell_circuits_scan/Base/Combo_CC_Base.py`),
`X = 0` and `Y = 0` are absorbing states and `X` is bounded by the carrying
capacity, so trajectories remain in the first quadrant and cannot diverge.

## Directory structure

```
Circuit_to_target_analysis/
  README.md
  circuit_to_target/                  # Main Python package
    model.py                          # ODEs, Euler integration, DynamicModel2D bridge
    metrics.py                        # threshold, pulse amplitude, refractory period
    sensitivity.py                    # logarithmic-derivative sensitivity and sweeps
    phase_diagram.py                  # dimensionless groups and regime boundaries
  notebooks/                          # Analysis pipeline (primary workflow)
    01_metric_panels.ipynb            # Fig. 6B to 6D
    02_sensitivity_analysis.ipynb     # Fig. 6E
    03_phase_diagram.ipynb            # Fig. 6F
  data/results/                       # Sensitivity tables written by notebook 02
  figures/                            # Figure panels written by the notebooks
```

## The three functional metrics

All three are defined in `circuit_to_target/metrics.py` and are integrated with
the explicit Euler method at `dt = 0.001`, using the non-negativity clamp of
`DynamicModel_Package.ModelBase.DynamicModel.euler_integrate_keep_positive`.

**Response threshold** (`threshold`) - the minimal `X` for which `dX/dt > 0` at
`Y = 0`. Solving the non-trivial `X` nullcline at `Y = 0`,
`(a c / d)(X/c)(1 - X/c) = 1`, gives the smaller root

```
X_threshold = c/2 - sqrt(a c (a c - 4 d)) / (2 a) = (c/2) [1 - sqrt(1 - 4 d / (c a))]
```

so the threshold depends only on `a`, `c` and `d`, and is exactly independent of
the inhibitor parameters `b`, `f` and `g`. When `a c < 4 d` the expression has no
real root: the `X` nullcline never reaches the first quadrant and no input can
excite a pulse.

**Pulse amplitude** (`pulse_amplitude`) - the maximal `X` reached from a
just-suprathreshold initial condition, with `Y(0) = 0.1`. Integration runs to
`t = 1` and, while the maximum still sits on the last time point, the window is
extended by a factor 1.5 so the true peak is always captured. The initial
effector level is set by `x0_rule`:

* `'threshold+1'` (default) - `X(0) = X_threshold + 1`;
* `'1.1*threshold'` - `X(0) = 1.1 X_threshold`.

The choice matters: at `1.1 X_threshold` the trajectory stays close to the
nullcline and the peak is set by the distance to threshold rather than by the
carrying capacity, which suppresses the dependence on `c` (see the two tables
below).

**Refractory period** (`refractory_period`) - the time during which `Y` stays
above 1% of its own maximum, integrated from `(X, Y) = (10, 0.1)` to `t = 10`,
extending the window while `Y` has not yet decayed. The maximum of `Y` is taken
over the whole integration.

## Sensitivity

`circuit_to_target/sensitivity.py` computes

```
S = (p / M) dM/dp = d log M / d log p
```

the percentage change in metric `M` per 1% change in parameter `p`. Each
parameter is varied individually over a logarithmically spaced range of 20
values while the others stay at the baseline set of Fig. 6E
(`a = b = f = g = 1`, `d = 2`, `c = 100`), and `S` is estimated by first-order
finite differences of the logarithms, averaged over the range.

Sweep ranges (`SWEEP_RANGES`):

| Metric | `a` | `b` | `c` | `d` | `f` | `g` |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| threshold | 0.1 - 10 | - | 20 - 500 | 0.01 - 20 | - | - |
| pulse amplitude | 0.5 - 2.5 | 0.1 - 5 | 20 - 500 | 0.01 - 10 | 0.1 - 2.5 | 0.1 - 10 |
| refractory period | 0.5 - 2.5 | 0.1 - 5 | 20 - 500 | 0.01 - 10 | 2.5 - 10 | 0.1 - 2 |

Parameters left blank do not enter the metric and are reported as exactly `0`.
For the threshold this is the analytic statement above; `log_sensitivity` also
requires the parameter and metric arrays to have equal length, so a sweep can
never be paired with the values of a different sweep.

Sensitivities produced by `02_sensitivity_analysis.ipynb`
(`data/results/sensitivity_table.csv`), with the default `x0_rule`:

| Metric | `a` | `b` | `c` | `d` | `f` | `g` |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| threshold | -1.070 | 0 | -0.036 | +1.043 | 0 | 0 |
| pulse amplitude | +1.734 | -0.813 | +1.106 | -0.185 | -0.968 | +0.246 |
| refractory period | -0.036 | -0.001 | -0.028 | +0.009 | -0.012 | -0.971 |

With `x0_rule='1.1*threshold'` the amplitude row instead reads
`a +1.892`, `b -0.905`, `c +0.008`, `d +1.101`, `f -1.170`, `g +0.813`
(`data/results/sensitivity_table_1.1threshold.csv`).

The threshold rises with the removal rate `d` and falls with the autocrine
strength `a`; the amplitude rises with `a` and with the carrying capacity `c`;
and the refractory period is controlled almost entirely by the inhibitor
turnover `g`, to which the other five parameters contribute little.

### Relation to the published Fig. 6E

The threshold and refractory-period panels are reproduced by this code up to
small numerical differences (largest bar heights agree to within 0.05; the
ordering, signs and conclusions are unchanged). Two differences are worth
recording:

* The published threshold panel used `d = 1` as the baseline for the `a` and `c`
  sweeps, rather than the `d = 2` of the Fig. 6E caption used here
  (`a`: -1.026 against -1.070).
* The published refractory panel compared `Y` against 1% of its maximum over the
  first `t = 1` of the integration, before `Y` peaks, rather than over the whole
  integration; this affects only the near-zero `c` bar (+0.066 against -0.028).

In the published pulse-amplitude panel the `b`, `f` and `g` bars were evaluated
against a placeholder parameter array rather than against the sweeps listed
above, and the `a` and `d` bars against the wider threshold-sweep ranges. Those
bars are therefore not reproduced by the table above: the amplitude of this
circuit does depend on the inhibitor parameters `b` and `f`, as the entries
`-0.813` and `-0.968` show.

## Dynamical regimes

`circuit_to_target/phase_diagram.py` non-dimensionalises the circuit with
`x~ = (a/d) X`, `y~ = (b/d) Y`, `t~ = d t`:

```
dx~/dt~ = x~ [ x~ (1 - x~/p1) - y~ - 1 ]
dy~/dt~ = y~ [ p2 x~ - p3 ]
```

with `p1 = c a / d` (self-activation strength), `p2 = f / a` (inhibitor
induction strength) and `p3 = g / d = 1` (no timescale separation). The
non-trivial `x~` nullcline meets the `x~` axis at
`X± = [p1 ± sqrt(p1 (p1 - 4))] / 2`, real only for `p1 >= 4`. This gives three
regimes:

| Regime | Condition | Interpretation |
| :--- | :--- | :--- |
| Monostable "OFF" | `p1 < 4` | no input excites a pulse; hypoimmunity |
| Bistable "ON"/"OFF" | `p1 > 4` and `p2 < 1 / X+` | pulse does not terminate; chronic activation, hyperimmunity |
| Excitable | `p1 > 4` and `p2 > 1 / X+` | suprathreshold input gives a pulse that returns to rest |

The middle boundary follows from the stability of the fixed point at `(X+, 0)`:
it is stable along `x~`, and along `y~` its Jacobian entry is `p2 X+ - p3`, so it
is a second attractor exactly when `p2 X+ < p3 = 1`. The boundary drawn in
Fig. 6F is therefore

```
p2 = 1 / X+ = 2 / [ p1 + sqrt(p1 (p1 - 4)) ]
```

`03_phase_diagram.ipynb` checks this classification against direct integration
of the full circuit. The circuit shows no oscillatory regime.

## Reproducibility checklist

The package needs `libs/` on the Python path only for `to_dynamic_model()`; see
the repository README:

```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)/../libs
```

1. **`01_metric_panels.ipynb`** - regenerates the illustrative panels of
   Fig. 6B to 6D (pulse amplitude for `a = 1` against `a = 5`; threshold for
   `a = 0.1` against `a = 3`; refractory period for `g = 1` against `g = 10`)
   and checks `integrate()` against `DynamicModel2D`.
2. **`02_sensitivity_analysis.ipynb`** - runs all sweeps, writes
   `data/results/sensitivity_table.csv` and the three Fig. 6E bar charts.
3. **`03_phase_diagram.ipynb`** - draws the Fig. 6F phase diagram, verifies the
   regime boundaries against direct integration, and saves one representative
   trajectory per regime.

## Dependencies

`numpy`, `pandas`, `matplotlib`, `jupyter`; `DynamicModel_Package` (bundled in
`libs/`) for `to_dynamic_model()`.
