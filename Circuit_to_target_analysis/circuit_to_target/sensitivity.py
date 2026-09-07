"""
Sensitivity analysis of the three circuit metrics
(Methods, "Circuit to target and sensitivity analysis"; Fig. 6E).

Sensitivity is the logarithmic derivative of a metric M with respect to a
parameter p,

    S = (p / M) dM/dp = d log M / d log p,

i.e. the percentage change in the metric per 1% change in the parameter. Each
parameter is varied individually over a logarithmically spaced range while the
others are held at the baseline set of Fig. 6E, and S is estimated by
first-order finite differences of the logarithms, averaged over the range.
"""

import numpy as np
import pandas as pd

from .metrics import pulse_amplitude, refractory_period, threshold
from .model import PARAMETERS_FIG6

#: Parameters of the circuit, in the order used by the Fig. 6E bar charts.
PARAMETER_NAMES = ('a', 'b', 'c', 'd', 'f', 'g')

#: Colours used per parameter in Fig. 6E.
PARAMETER_COLORS = ('tab:blue', 'tab:orange', 'tab:green',
                    'tab:pink', 'tab:red', 'tab:brown')

#: Range over which each parameter is swept, per metric. The threshold is
#: analytic and is swept over a wider range than the two dynamical metrics.
SWEEP_RANGES = {
    'threshold':       {'a': (0.1, 10.0), 'c': (20.0, 500.0), 'd': (0.01, 20.0)},
    'pulse_amplitude': {'a': (0.5, 2.5), 'b': (0.1, 5.0), 'c': (20.0, 500.0),
                        'd': (0.01, 10.0), 'f': (0.1, 2.5), 'g': (0.1, 10.0)},
    'refractory':      {'a': (0.5, 2.5), 'b': (0.1, 5.0), 'c': (20.0, 500.0),
                        'd': (0.01, 10.0), 'f': (2.5, 10.0), 'g': (0.1, 2.0)},
}

#: Number of parameter values per sweep.
N_SWEEP = 20


def log_sensitivity(parameter_values, metric_values):
    """
    Estimate d log M / d log p by first-order finite differences.

    `parameter_values` and `metric_values` must be the same length and must
    describe the *same* sweep: pairing a metric computed over one parameter
    range with a different range of parameter values silently returns a
    rescaled sensitivity, so the lengths are checked here.
    """
    p = np.asarray(parameter_values, dtype=float)
    m = np.asarray(metric_values, dtype=float)
    if p.shape != m.shape:
        raise ValueError(
            f"parameter_values has length {p.size} but metric_values has "
            f"length {m.size}; they must describe the same sweep")
    return float(np.mean(np.diff(np.log(m)) / np.diff(np.log(p))))


def sweep(parameter, low, high, n=N_SWEEP):
    """Logarithmically spaced values of one parameter."""
    return np.logspace(np.log10(low), np.log10(high), n)


def sweep_metric(metric, parameter, low, high, baseline=None, n=N_SWEEP,
                 **metric_kwargs):
    """
    Evaluate `metric` over a sweep of one parameter, holding the others at
    `baseline` (default: the Fig. 6E parameter set).

    `metric` is called as ``metric(parameters, **metric_kwargs)``, except for
    `threshold`, which takes a, c and d directly.

    Returns
    -------
    (values, metric_values) : tuple of np.ndarray
    """
    baseline = dict(PARAMETERS_FIG6 if baseline is None else baseline)
    values = sweep(parameter, low, high, n)

    out = []
    for value in values:
        pars = dict(baseline)
        pars[parameter] = value
        if metric is threshold:
            out.append(threshold(pars['a'], pars['c'], pars['d']))
        else:
            out.append(metric(pars, **metric_kwargs))
    return values, np.asarray(out, dtype=float)


def sensitivity_table(baseline=None, n=N_SWEEP, x0_rule='threshold+1'):
    """
    Sensitivity of all three metrics to all six parameters (Fig. 6E).

    Parameters not listed in `SWEEP_RANGES` for a given metric do not affect
    it and are reported as exactly 0. For the threshold this is an analytic
    statement: b, f and g do not appear in the threshold expression at all.

    Returns
    -------
    pd.DataFrame
        Rows 'threshold', 'pulse_amplitude', 'refractory'; columns a...g.
    """
    metrics = {
        'threshold': (threshold, {}),
        'pulse_amplitude': (pulse_amplitude, {'x0_rule': x0_rule}),
        'refractory': (refractory_period, {}),
    }

    rows = {}
    for name, (metric, kwargs) in metrics.items():
        row = {}
        for parameter in PARAMETER_NAMES:
            if parameter not in SWEEP_RANGES[name]:
                row[parameter] = 0.0
                continue
            low, high = SWEEP_RANGES[name][parameter]
            values, metric_values = sweep_metric(
                metric, parameter, low, high, baseline=baseline, n=n, **kwargs)
            row[parameter] = log_sensitivity(values, metric_values)
        rows[name] = row

    return pd.DataFrame(rows).T[list(PARAMETER_NAMES)]
