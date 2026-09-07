"""
The three functional metrics of the circuit-to-target analysis
(Methods, "Circuit to target and sensitivity analysis"; Fig. 6B to 6D):

    * `threshold`         - response threshold: minimal X that excites a pulse
    * `pulse_amplitude`   - response amplitude: maximal X reached in the pulse
    * `refractory_period` - time until the inhibitor Y falls back to 1% of its
                            maximum

Each metric is a scalar function of the six circuit parameters, so that
`sensitivity.py` can take its logarithmic derivative with respect to each
parameter in turn.
"""

import numpy as np

from .model import DT, integrate

#: Inhibitor level the pulse-amplitude integration starts from.
Y0_AMPLITUDE = 0.1

#: Initial condition for the refractory-period integration (Methods: X = 10,
#: Y = 0.1).
X0_REFRACTORY = 10.0
Y0_REFRACTORY = 0.1

#: Fraction of the inhibitor peak that marks the end of the refractory period.
REFRACTORY_FRACTION = 0.01


def threshold(a, c, d, **ignored):
    """
    Response threshold: the smallest X for which dX/dt > 0 at Y = 0.

    Setting the non-trivial X nullcline to zero,

        (a c / d) (X / c) (1 - X / c) = 1,

    whose smaller root is

        X_threshold = c/2 - sqrt(a c (a c - 4 d)) / (2 a)
                    = (c/2) [1 - sqrt(1 - 4 d / (c a))].

    The threshold therefore depends only on a, c and d; it is exactly
    independent of the inhibitor parameters b, f and g, which is why those
    three sensitivities are analytically zero (see `sensitivity.py`).

    Returns NaN when a c < 4 d, i.e. when the X nullcline never reaches the
    first quadrant and no pulse is possible ("OFF" state).
    """
    discriminant = 1.0 - 4.0 * d / (c * a)
    if discriminant < 0:
        return np.nan
    return 0.5 * c * (1.0 - np.sqrt(discriminant))


def pulse_amplitude(parameters, x0_rule='threshold+1', y0=Y0_AMPLITUDE,
                    t_final=1.0, dt=DT, max_extensions=40):
    """
    Response amplitude: the maximal X reached from a just-suprathreshold start.

    The circuit is integrated forward from an initial condition slightly above
    the threshold for the given parameter set. Integration runs to `t_final`
    and, if X is still rising at the end of the window (the maximum sits on the
    last time point), the window is extended by a factor 1.5 and the
    integration repeated, so that the true peak is always captured.

    `x0_rule` selects the initial effector level:

    ``'threshold+1'``
        X(0) = X_threshold + 1. This is the rule used for Fig. 6E.
    ``'1.1*threshold'``
        X(0) = 1.1 X_threshold, the "110% threshold" wording of the Methods.

    The two rules do not give the same sensitivities - at 1.1 X_threshold the
    trajectory stays close to the nullcline and the peak is set by the distance
    to threshold rather than by the carrying capacity, which suppresses the
    dependence on c. See the repository README for the comparison.
    """
    x_threshold = threshold(parameters['a'], parameters['c'], parameters['d'])
    if not np.isfinite(x_threshold):
        return np.nan

    if x0_rule == 'threshold+1':
        x0 = x_threshold + 1.0
    elif x0_rule == '1.1*threshold':
        x0 = 1.1 * x_threshold
    else:
        raise ValueError(f"unknown x0_rule {x0_rule!r}; expected "
                         "'threshold+1' or '1.1*threshold'")

    t = t_final
    for _ in range(max_extensions):
        traj = integrate(x0, y0, parameters, t, dt=dt)
        if traj['x'].values.argmax() != len(traj) - 1:
            break
        t *= 1.5
    return float(traj['x'].max())


def refractory_period(parameters, x0=X0_REFRACTORY, y0=Y0_REFRACTORY,
                      fraction=REFRACTORY_FRACTION, t_final=10.0, dt=DT,
                      max_extensions=40):
    """
    Refractory period: the time during which the inhibitor Y stays above
    `fraction` of its own maximum.

    The circuit is integrated from (X, Y) = (10, 0.1). The peak of Y is found,
    and the refractory period is the last time at which Y is still above
    `fraction` x max(Y) - equivalently the first time after the peak at which Y
    drops below it. If Y has not yet decayed by the end of the window, the
    window is extended by a factor 1.5 and the integration repeated.
    """
    t = t_final
    for _ in range(max_extensions):
        traj = integrate(x0, y0, parameters, t, dt=dt)
        ys = traj['y'].values
        if ys[-1] < ys.max() * fraction:
            break
        t *= 1.5

    ys = traj['y'].values
    above = np.nonzero(ys > ys.max() * fraction)[0]
    if len(above) == 0:
        return 0.0
    return float(traj['t'].values[above[-1]])
