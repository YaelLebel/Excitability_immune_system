"""
Dynamical regimes of the reduced excitable circuit
(Methods, "Circuit to target and sensitivity analysis" - Phase diagram; Fig. 6F).

Non-dimensionalising the circuit of `model.py` with

    x~ = (a/d) X,   y~ = (b/d) Y,   t~ = d t

removes three of the six parameters and leaves

    dx~/dt~ = x~ [ x~ (1 - x~/p1) - y~ - 1 ]
    dy~/dt~ = y~ [ p2 x~ - p3 ]

with the dimensionless groups

    p1 = c a / d    self-activation strength
    p2 = f / a      inhibitor induction strength
    p3 = g / d      inhibitor turnover, set to 1 (no timescale separation).

Regimes
-------
The non-trivial x~ nullcline is y~ = x~ (1 - x~/p1) - 1, an inverted parabola
that meets the x~ axis where x~^2 - p1 x~ + p1 = 0, i.e. at

    X± = [ p1 ± sqrt(p1 (p1 - 4)) ] / 2,

which are real only for p1 >= 4.

* ``p1 < 4``  - the nullcline never reaches the first quadrant, so no input can
  excite a pulse and the origin is the only stable fixed point: monostable
  "OFF" (hypoimmunity).

* ``p1 > 4`` - the fixed point at (X+, 0) is stable along x~ (the parabola is
  positive between the roots and negative beyond X+). Its stability along y~ is
  set by the sign of dy~/dt~ there, whose Jacobian entry is p2 X+ - p3. With
  p3 = 1 it is stable exactly when

      p2 X+ < 1,   i.e.   p2 < 1 / X+ = 2 / [ p1 + sqrt(p1 (p1 - 4)) ].

  In that case (X+, 0) is a second stable fixed point and the circuit is
  bistable "ON"/"OFF" (chronic activation, hyperimmunity). Above that boundary
  the inhibitor is induced strongly enough to destabilise (X+, 0), the origin is
  the only attractor, and a suprathreshold input produces a pulse that returns
  to rest: excitable.

So the excitable regime is p1 > 4 together with
p2 > 2 / [p1 + sqrt(p1 (p1 - 4))], and the boundary drawn in Fig. 6F is
p2 = 1 / X+.
"""

import numpy as np
from matplotlib.ticker import (FixedLocator, FuncFormatter,
                               NullLocator)

#: Value of p1 below which the x~ nullcline does not reach the first quadrant.
P1_CRITICAL = 4.0

REGIME_COLORS = {
    'monostable_off': 'tab:blue',
    'bistable': 'tab:green',
    'excitable': 'tab:orange',
}


def dimensionless_groups(parameters):
    """
    Map the six circuit parameters onto (p1, p2, p3).

    p1 = c a / d, p2 = f / a, p3 = g / d.
    """
    a, c, d = parameters['a'], parameters['c'], parameters['d']
    return {
        'p1': c * a / d,
        'p2': parameters['f'] / a,
        'p3': parameters['g'] / d,
    }


def x_plus(p1):
    """
    Larger root of the x~ nullcline, X+ = [p1 + sqrt(p1 (p1 - 4))] / 2.

    Returns NaN for p1 < 4, where the roots are complex.
    """
    p1 = np.asarray(p1, dtype=float)
    with np.errstate(invalid='ignore'):
        root = np.sqrt(p1 * (p1 - P1_CRITICAL))
    return np.where(p1 >= P1_CRITICAL, 0.5 * (p1 + root), np.nan)


def bistable_boundary(p1):
    """
    The p2 below which (X+, 0) is stable and the circuit is bistable:

        p2 = 1 / X+ = 2 / [p1 + sqrt(p1 (p1 - 4))].
    """
    return 1.0 / x_plus(p1)


def classify_regime(p1, p2, p3=1.0):
    """
    Return 'monostable_off', 'bistable' or 'excitable' for one (p1, p2) pair.

    `p3` is accepted for completeness; the boundaries above assume p3 = 1, as
    in Fig. 6F.
    """
    if p1 < P1_CRITICAL:
        return 'monostable_off'
    if p2 * x_plus(p1) < p3:
        return 'bistable'
    return 'excitable'


def plot_phase_diagram(ax, p1_lim=(2.0, 15.0), p2_lim=(0.1, 15.0), n=250,
                       alpha=0.5, lw=3):
    """
    Draw the (p1, p2) phase diagram of Fig. 6F on `ax`.

    The two boundaries are the vertical line p1 = 4 and the curve
    p2 = 1 / X+(p1); the three regions are shaded with `REGIME_COLORS`.
    """
    p1_min, p1_max = p1_lim
    p2_min, p2_max = p2_lim

    # Monostable "OFF": p1 < 4, all p2.
    ax.fill_betweenx([p2_min, p2_max], [p1_min] * 2, [P1_CRITICAL] * 2,
                     alpha=alpha, color=REGIME_COLORS['monostable_off'])
    ax.plot([P1_CRITICAL] * 2, [p2_min, p2_max], lw=lw, c='black')

    # Boundary between bistable and excitable, for p1 > 4.
    p1s = np.linspace(P1_CRITICAL * 1.0025, p1_max, n)
    boundary = bistable_boundary(p1s)
    ax.plot(p1s, boundary, lw=lw, c='black')

    ax.fill_between(p1s, p2_min, boundary, alpha=alpha,
                    color=REGIME_COLORS['bistable'])
    ax.fill_between(p1s, boundary, p2_max, alpha=alpha,
                    color=REGIME_COLORS['excitable'])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([p1_min, p1_max])
    ax.set_ylim([p2_min, p2_max])

    # Plain tick labels: the default log minor labels collide on this range.
    ax.xaxis.set_major_locator(FixedLocator([2, 4, 6, 10]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    ax.yaxis.set_major_locator(FixedLocator([0.1, 1, 10]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    ax.set_xlabel('$p_1$ (self-activation strength)', fontsize=12)
    ax.set_ylabel('$p_2$ (inhibitor induction strength)', fontsize=12)
    return ax
