"""
Numerical integration utilities for `Combo` circuits (polynomial ODEs
without explicit carrying capacities).

`Combo_ODE` wraps a symbolic `Combo` object and:
- builds a right‑hand side function suitable for `scipy.integrate.odeint`,
- substitutes specific parameter values (dp0, dp1, …),
- returns trajectories as a pandas DataFrame,
- optionally extends the integration time until both X and Y have relaxed
  back after an excitable pulse (via `continue_to_decay`).
"""

import numpy as np
import sympy as sp
from scipy.integrate import odeint
from utils import cartes_to_polar
import matplotlib.pyplot as plt
import pandas as pd
class Combo_ODE ():
    def __init__ (self,combo,parameters_values,tol = 1e-4):
        self.combo = combo
        self.parameters_values = parameters_values
        self.tol = tol
    
    def fun_to_integrate (self,Y,t):
        x,y = Y
        dxdt = self.combo.P
        if len(self.parameters_values)>0:
            for p in self.combo.params:
                dxdt = dxdt.subs(p,self.parameters_values[p])
        dxdt = dxdt.subs({self.combo.x:x,self.combo.y:y})
        dydt = self.combo.Q
        if len(self.parameters_values)>0:
            for p in self.combo.params:
                dydt = dydt.subs(p,self.parameters_values[p])
        dydt = dydt.subs({self.combo.x:x,self.combo.y:y})
        return [dxdt,dydt]
    

    def integrate (self,Y0,t_final,dt,continue_to_decay = False):
        t_final= t_final
        dt = dt
        t = np.arange(0,t_final+dt,dt)
        ret = odeint(self.fun_to_integrate,Y0,t)
        xs = ret[:,0]
        ys = ret[:,1]
        if continue_to_decay:
            flag_cont = False
            while not flag_cont:
                if (np.argmax(xs)==len(xs)-1 and xs[-1]-xs[-2]>self.tol) or (np.argmax(ys)==len(ys)-1 and ys[-1]-ys[-2]>self.tol):
                    t_final = 1.5*t_final
                    ret = self.integrate(Y0,t_final,dt)
                    xs = ret[:,0]
                    ys = ret[:,1]
                else:
                    flag_cont = True
        df_ret = pd.DataFrame({'t':t,'x':xs,'y':ys})
        return df_ret
    
    def integrate_and_plot (self,Y0,t_final,dt,fig = None, axs = None,continue_to_decay = False):
        ret = self.integrate(Y0,t_final,dt)
        xs = ret[:,0]
        ys = ret[:,1]
        if continue_to_decay:
            flag_cont = False
            while not flag_cont:
                if (np.argmax(xs)==len(xs)-1) or (np.argmax(ys)==len(ys)-1):
                    t_final = 1.5*t_final
                    ret = self.integrate(Y0,t_final,dt)
                    xs = ret[:,0]
                    ys = ret[:,1]
                else:
                    flag_cont = True
        flag_return_fig = False
        if fig == None and axs == None:
            fig,axs = plt.subplots(nrows = 2, ncols = 1)
            flag_return_fig = True
        axs[0].plot(self.t,ret[:,0])
        axs[0].set_title('x')
        axs[1].plot(self.t,ret[:,1])
        axs[1].set_title('y')
        fig.tight_layout()
        if flag_return_fig:
            return ret, fig, axs
        else:
            return ret


