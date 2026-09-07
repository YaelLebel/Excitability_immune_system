import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from .circuit import Circuit


class CircuitPlotter:
    """
    Handles plotting of phase portraits and nullclines.
    Refactored from Combo_plot.

    The vector field and the nullclines are evaluated numerically with
    `sympy.lambdify`, so this class is self-contained and does not depend on
    `DynamicModel_Package.DynamicModel2D`.
    """

    def __init__(self, circuit: Circuit, parameters_dict: dict = {}):
        self.circuit = circuit
        self.parameters_dict = parameters_dict
        self.x, self.y = self.circuit.variables
        self.parameter_values = list(self.parameters_dict.values())

    def deriv_x(self, vars: dict, pars: dict):
        """Evaluate dx/dt at the point given by `vars`, with parameters `pars`."""
        dict_all = vars.copy()
        if len(pars) > 0:
            dict_all.update(pars)
        return self.circuit.xdot.subs(dict_all)

    def deriv_y(self, vars: dict, pars: dict):
        """Evaluate dy/dt at the point given by `vars`, with parameters `pars`."""
        dict_all = vars.copy()
        if len(pars) > 0:
            dict_all.update(pars)
        return self.circuit.ydot.subs(dict_all)

    def det_y_fun_x(self, n):
        """Whether nullcline branch `n` is stored as y(x) (True) or x(y) (False)."""
        y_fun_x = True
        if str(n[1]) == 'y':
            y_fun_x = False
        return y_fun_x

    def plot_phase_portrait(self, t, n_X, n_Y, xlim, ylim, colors=None,
                            fig=None, ax=None, density=1.0, field_color='grey'):
        flag_return_figs = False
        if fig is None and ax is None:
            fig, ax = plt.subplots()
            flag_return_figs = True

        # Nullclines are optional in Circuit.__init__; solve them on demand so
        # the portrait is never silently drawn without them.
        if not self.circuit.nullclines:
            self.circuit.find_nullclines()

        # Vector field
        x_vals = np.linspace(xlim[0], xlim[1], n_X)
        y_vals = np.linspace(ylim[0], ylim[1], n_Y)
        X, Y = np.meshgrid(x_vals, y_vals)

        dxdt_expr = self.circuit.xdot.subs(self.parameters_dict)
        dydt_expr = self.circuit.ydot.subs(self.parameters_dict)

        f_x = sp.lambdify((self.x, self.y), dxdt_expr, modules='numpy')
        f_y = sp.lambdify((self.x, self.y), dydt_expr, modules='numpy')

        # lambdify collapses a component with no x/y dependence to a scalar,
        # so broadcast back to the grid shape before streamplot.
        U = np.broadcast_to(np.asarray(f_x(X, Y), dtype=float), X.shape).copy()
        V = np.broadcast_to(np.asarray(f_y(X, Y), dtype=float), X.shape).copy()

        # Grey, behind the nullclines, as in DynamicModel2D.plot_streamplot:
        # the default streamline colour clashes with the x nullcline.
        ax.streamplot(X, Y, U, V, density=density, color=field_color, zorder=0)

        # Nullclines
        self._plot_nullclines(self.circuit.nullclines_x, xlim, ylim,
                              color=colors[0] if colors else 'tab:blue',
                              label='x nullcline', ax=ax)
        self._plot_nullclines(self.circuit.nullclines_y, xlim, ylim,
                              color=colors[1] if colors else 'tab:orange',
                              label='y nullcline', ax=ax)

        ax.legend(loc='upper right')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(str(self.x))
        ax.set_ylabel(str(self.y))

        if flag_return_figs:
            return fig, ax
        return fig, ax

    def _nullcline_points(self, nc, xlim, ylim, n_points=200):
        """
        Evaluate one nullcline branch on a grid, keeping only the points that
        are real, finite and positive in the dependent variable.

        Returns (xs, ys) ready to plot.
        """
        y_fun_x = self.det_y_fun_x(nc)
        expr = (nc[1] if y_fun_x else nc[0]).subs(self.parameters_dict)
        free_var = self.x if y_fun_x else self.y
        f = sp.lambdify(free_var, expr, modules='numpy')

        lim = xlim if y_fun_x else ylim
        indep = np.linspace(lim[0], lim[1], n_points)

        dep = np.broadcast_to(np.asarray(f(indep), dtype=complex), indep.shape)

        # Discard complex and non-finite branches before comparing to 0:
        # NumPy raises on order comparisons of complex values.
        keep = np.isclose(dep.imag, 0.0) & np.isfinite(dep.real)
        indep, dep = indep[keep], dep[keep].real

        keep = dep > 0
        indep, dep = indep[keep], dep[keep]

        return (indep, dep) if y_fun_x else (dep, indep)

    def _plot_nullclines(self, nullclines, xlim, ylim, color, label, ax):
        """Plot every branch of one nullcline, labelling the legend entry once."""
        labelled = False
        for nc in nullclines:
            xs, ys = self._nullcline_points(nc, xlim, ylim)
            if len(xs) == 0:
                continue
            ax.plot(xs, ys, lw=3, color=color, zorder=1,
                    label=None if labelled else label)
            labelled = True
