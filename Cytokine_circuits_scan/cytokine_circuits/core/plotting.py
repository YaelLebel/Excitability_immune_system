import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from .circuit import Circuit
# from DynamicModel_Package.DynamicModel2D_Base import DynamicModel2D 
# Note: DynamicModel2D_Base was used in original but source not visible in Base dir.
# It was imported from 'DynamicModel_Package'.
# I'll check if I need to implement streamplot logic manually if that package isn't available.

class CircuitPlotter:
    """
    Handles plotting of phase portraits and nullclines.
    Refactored from Combo_plot.
    """
    def __init__(self, circuit: Circuit, parameters_dict: dict = {}):
        self.circuit = circuit
        self.parameters_dict = parameters_dict
        self.x, self.y = self.circuit.variables
        self.model = self.create_model()
        self.parameter_values = list(self.parameters_dict.values())

    
    def deriv_x (self, vars:dict,pars:dict):
        if len(pars)>0:
            dict_all = vars.copy()
            dict_all.update(pars)
        else:
            dict_all = vars
        expr_str = str(self.combo.xdot)
        expr = sp.sympify(expr_str)
        ret = expr.subs(dict_all)
        return ret
    
    def deriv_y (self, vars:dict,pars:dict):
        if len(pars)>0:
            dict_all = vars.copy()
            dict_all.update(pars)
        else:
            dict_all = vars
        expr_str = str(self.combo.ydot)
        expr = sp.sympify(expr_str)
        ret = expr.subs(dict_all)
        return ret
    
    def create_model(self):
        model = DynamicModel2D('x',self.deriv_x,self.parameters_dict,
                               'y',self.deriv_y,self.parameters_dict)
        return model

    def det_y_fun_x (self,n):
        y_fun_x = True
        if str(n[1])=='y':
            y_fun_x = False
        return y_fun_x

    def plot_phase_portrait(self, t, n_X, n_Y, xlim, ylim, colors=None, fig=None, ax=None):

        if fig is None and ax is None:
            fig, ax = plt.subplots()
            
        # Streamplot logic (replaces DynamicModel2D dependency if needed)
        x_vals = np.linspace(xlim[0], xlim[1], n_X)
        y_vals = np.linspace(ylim[0], ylim[1], n_Y)
        X, Y = np.meshgrid(x_vals, y_vals)
        
        # Calculate vector field
        # Lambdify for speed
        dxdt_expr = self.circuit.xdot.subs(self.parameters_dict)
        dydt_expr = self.circuit.ydot.subs(self.parameters_dict)
        
        f_x = sp.lambdify((self.x, self.y), dxdt_expr, modules='numpy')
        f_y = sp.lambdify((self.x, self.y), dydt_expr, modules='numpy')
        
        U = f_x(X, Y)
        V = f_y(X, Y)
        
        ax.streamplot(X, Y, U, V, density=1.0)
        
        # Plot Nullclines
        self._plot_nullclines(self.circuit.nullclines_x, xlim, ylim, color=(colors[0] if colors else 'tab:blue'), label='x nullcline', ax=ax)
        self._plot_nullclines(self.circuit.nullclines_y, xlim, ylim, color=(colors[1] if colors else 'tab:orange'), label='y nullcline', ax=ax)

        ax.legend(loc='upper right')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        
        return fig, ax

    def _plot_nullclines(self, nullclines, xlim, ylim, color, label, ax):
        # Helper to plot nullclines
        for nc in nullclines:            
            if self.det_y_fun_x(nc):
                nc_f = sp.lambdify(x,nc[1].subs(self.parameters_dict))
                x_vs_nc = np.linspace(xlim[0],xlim[1],200)
                y_vs_nc = nc_f(x_vs_nc)
                mask_positive = y_vs_nc>0
                x_vs_nc = x_vs_nc[mask_positive]
                y_vs_nc = y_vs_nc[mask_positive]
                mask_real = np.imag(y_vs_nc)==0
                x_vs_nc = x_vs_nc[mask_real]
                y_vs_nc = y_vs_nc[mask_real]
                if colors == None:
                    ax.plot(x_vs_nc,y_vs_nc,lw = 3,label = 'x nullcline')
                else:
                    ax.plot(x_vs_nc,y_vs_nc,lw = 3,label = 'x nullcline',color = colors[0])
            else:
                nc_f = sp.lambdify(y,nc[0].subs(self.parameters_dict))
                y_vs_nc = np.linspace(ylim[0],ylim[1],200)
                x_vs_nc = nc_f(y_vs_nc)
                mask_positive = x_vs_nc>0
                x_vs_nc = x_vs_nc[mask_positive]
                y_vs_nc = y_vs_nc[mask_positive]
                mask_real = np.imag(x_vs_nc)==0
                x_vs_nc = x_vs_nc[mask_real]
                y_vs_nc = y_vs_nc[mask_real]

                if colors == None:
                    ax.plot(x_vs_nc,y_vs_nc,lw = 3,label = 'x nullcline')
                else:
                    ax.plot(x_vs_nc,y_vs_nc,lw = 3,label = 'x nullcline',color = colors[0])
