import numpy as np
import sympy as sp
from Base.Combo_pp_Base import Combo_pp
from DynamicModel_Package.DynamicModel2D_Base import DynamicModel2D
import sympy as sp
from sympy.utilities.lambdify import lambdastr
import matplotlib.pyplot as plt
x,y = sp.symbols('x y')

class Combo_plot ():
    def __init__ (self,combo:Combo_pp,parameters_dict:dict = {}):
        self.combo = combo
        self.parameters_dict = parameters_dict
        self.model = self.create_model()
        self.parameter_values = list(self.parameters_dict.values())
        #self.convert_nullclines()

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


    def plot_phase_portrait(self,t:float,n_X:int,n_Y:int, xlim:tuple, ylim:tuple,colors:list = None, fig = None, ax = None):
        flag_return_figs = False
        if fig == None and ax == None:
            fig,ax = self.model.plot_streamplot(t,n_X,n_Y,xlim,ylim,fig,ax)
            flag_return_figs = True
        else:
            self.model.plot_streamplot(t,n_X,n_Y,xlim,ylim,fig,ax)
        for nc in self.combo.nullclines_x:
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
        for nc in self.combo.nullclines_y:
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
                    ax.plot(x_vs_nc,y_vs_nc,lw = 3,label = 'y nullcline')
                else:
                    ax.plot(x_vs_nc,y_vs_nc,lw = 3,label = 'y nullcline',color = colors[1])
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
                    ax.plot(x_vs_nc,y_vs_nc,lw = 3,label = 'y nullcline')
                else:
                    ax.plot(x_vs_nc,y_vs_nc,lw = 3,label = 'y nullcline',color = colors[1])
        ax.legend(loc = 'upper right')
        ax.set_xlim(xlim[0],xlim[1])
        ax.set_ylim(ylim[0],ylim[1])
        if flag_return_figs:
            return fig,ax
        

"""


    def convert_nullcline(self,expr,variable):
        expr_str = str(expr)
        expr = sp.sympify(expr_str)

        def f (var,pars):
            if len(pars)>0:
                if str(variable)=='x':
                    new_dict = {'x':var}
                    new_dict.update(pars)
                else:
                    new_dict = {'y':var}
                    new_dict.update(pars)
            else:
                if str(variable)=='x':
                    new_dict = {'x':var}
                else:
                    new_dict = {'y':var}
            return expr.subs(new_dict)

        return f

    def convert_nullclines(self):
        for n in self.combo.nullclines_x:
            y_fun_x = self.det_y_fun_x(n)
            if y_fun_x:
                f = self.convert_nullcline(n[1],self.combo.x)
            else:
                f = self.convert_nullcline(n[0],self.combo.y)
            self.model.add_nullcline('x',f,self.parameters_dict,y_fun_x)
        for n in self.combo.nullclines_y:
            y_fun_x = self.det_y_fun_x(n)
            if y_fun_x:
                f = self.convert_nullcline(n[1],self.combo.x)
            else:
                f = self.convert_nullcline(n[0],self.combo.y)
            self.model.add_nullcline('y',f,self.parameters_dict,y_fun_x)

    def plot_phase_portrait(self,t:float,n_X:int,n_Y:int, xlim:tuple, ylim:tuple,colors:list = None, fig = None, ax = None):
        return self.model.plot_phase_portrait(t,n_X,n_Y,xlim,ylim,colors,fig,ax)
    


"""