"""
Reduced excitable cell-cell circuit used for the circuit-to-target analysis
(Methods, "Circuit to target and sensitivity analysis"; Fig. 6).

This is the reduced form of the excitable architecture identified by the
five-interaction cell-cell screen in `Cell_circuits_scan/`:

    dX/dt = X [ a X (1 - X/c) - b Y - d ]
    dY/dt = Y [ f X - g ]

with the six parameters of Fig. 6A:

    a   autocrine (self-activation) strength of the effector X   [1 / (X t)]
    b   inhibition of X by the inhibitor Y                       [1 / (Y t)]
    c   carrying capacity of X                                   [X]
    d   removal rate of X                                        [1 / t]
    f   induction of the inhibitor Y by X                        [1 / (X t)]
    g   turnover (removal) rate of Y                             [1 / t]

As in the general cell-cell class (`Cell_circuits_scan/Base/Combo_CC_Base.py`),
both X = 0 and Y = 0 are absorbing states and X is bounded by the carrying
capacity c, so trajectories stay in the first quadrant and cannot diverge.

Integration uses the explicit Euler scheme with the non-negativity clamp of
`DynamicModel_Package.ModelBase.DynamicModel.euler_integrate_keep_positive`
(any value below `tol` is set to 0). It is reimplemented here directly on NumPy
arrays because each sensitivity sweep integrates the circuit several hundred
times; `to_dynamic_model()` returns the equivalent `DynamicModel2D` object used
elsewhere in this repository, and
`notebooks/01_metric_panels.ipynb` checks that the two agree.
"""

import numpy as np
import pandas as pd

#: Integration time step used throughout Methods, "Circuit to target and
#: sensitivity analysis" (Euler, dt = 0.001).
DT = 1e-3

#: Values below this are clamped to 0, matching
#: `DynamicModel.euler_integrate_keep_positive`.
TOL = 1e-14

#: Baseline parameter set that the sensitivity analysis is taken around
#: (Fig. 6E caption).
PARAMETERS_FIG6 = {'a': 1.0, 'b': 1.0, 'c': 100.0, 'd': 2.0, 'f': 1.0, 'g': 1.0}


def dxdt(x, y, a, b, c, d, f, g):
    """Effector rate of change, dX/dt = X [a X (1 - X/c) - b Y - d]."""
    return x * (a * x * (1.0 - x / c) - b * y - d)


def dydt(x, y, a, b, c, d, f, g):
    """Inhibitor rate of change, dY/dt = Y [f X - g]."""
    return y * (f * x - g)


def integrate(x0, y0, parameters, t_final, dt=DT, tol=TOL):
    """
    Integrate the circuit forward in time with the explicit Euler method.

    Parameters
    ----------
    x0, y0 : float
        Initial effector and inhibitor levels.
    parameters : dict
        Values for 'a', 'b', 'c', 'd', 'f' and 'g'.
    t_final : float
        End of the integration window.
    dt : float, optional
        Time step (default 0.001, as in the Methods).
    tol : float, optional
        Values below this are set to 0, keeping both populations non-negative.

    Returns
    -------
    pd.DataFrame
        Columns 't', 'x' and 'y', one row per time step.
    """
    a, b, c, d, f, g = (parameters[k] for k in ('a', 'b', 'c', 'd', 'f', 'g'))
    n_steps = int(round(t_final / dt))

    xs = np.empty(n_steps + 1)
    ys = np.empty(n_steps + 1)
    x, y = float(x0), float(y0)
    xs[0], ys[0] = x, y

    for i in range(1, n_steps + 1):
        x_next = x + dxdt(x, y, a, b, c, d, f, g) * dt
        y_next = y + dydt(x, y, a, b, c, d, f, g) * dt
        x = 0.0 if x_next < tol else x_next
        y = 0.0 if y_next < tol else y_next
        xs[i], ys[i] = x, y

    return pd.DataFrame({'t': np.arange(n_steps + 1) * dt, 'x': xs, 'y': ys})


def to_dynamic_model(parameters):
    """
    Return this circuit as a `DynamicModel2D`, the shared two-variable model
    class used by `Cell_circuits_scan/` and `Cytokine_circuits_scan/`.

    Requires `libs/` on the Python path (see the repository README). Useful for
    phase portraits and for cross-checking `integrate()`.
    """
    from DynamicModel_Package.DynamicModel2D_Base import DynamicModel2D

    def _dx(variables, pars):
        return dxdt(variables['x'], variables['y'],
                    pars['a'], pars['b'], pars['c'],
                    pars['d'], pars['f'], pars['g'])

    def _dy(variables, pars):
        return dydt(variables['x'], variables['y'],
                    pars['a'], pars['b'], pars['c'],
                    pars['d'], pars['f'], pars['g'])

    return DynamicModel2D('x', _dx, parameters, 'y', _dy, parameters)
