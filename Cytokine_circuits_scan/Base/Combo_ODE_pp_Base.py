import numpy as np
import sympy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
from Base.Combo_pp_Base import Combo_pp
class Combo_ODE_pp ():
    def __init__ (self,combo:Combo_pp,parameters_values):
        self.combo = combo
        self.parameters_values = parameters_values
    
    def fun_to_integrate (self,Y,t):
        x,y = Y
        dxdt = self.combo.xdot
        if len(self.parameters_values)>0:
            for p in self.combo.parameters:
                dxdt = dxdt.subs(p,self.parameters_values[p])
        dxdt = dxdt.subs({self.combo.x:x,self.combo.y:y})
        dydt = self.combo.ydot
        if len(self.parameters_values)>0:
            for p in self.combo.parameters:
                dydt = dydt.subs(p,self.parameters_values[p])
        dydt = dydt.subs({self.combo.x:x,self.combo.y:y})
        return [dxdt,dydt]
    

    def integrate (self,Y0,t_final,dt,continue_to_decay = False,tol = 1e-4):
        t_final= t_final
        dt = dt
        t = np.arange(0,t_final+dt,dt)
        ret = odeint(self.fun_to_integrate,Y0,t)
        xs = ret[:,0]
        ys = ret[:,1]
        if continue_to_decay:
            flag_cont = False
            while not flag_cont:
                if abs(xs[-1]-xs[-2])>tol or abs(ys[-1]-ys[-2])>tol:
                    t_final = 1.5*t_final
                    ret = self.integrate(Y0,t_final,dt)
                    t = ret['t'].values
                    xs = ret['x'].values
                    ys = ret['y'].values
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


