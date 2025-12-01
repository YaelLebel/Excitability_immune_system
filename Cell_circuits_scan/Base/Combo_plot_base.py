import numpy as np
import sympy as sp
from coeff_base import coeff
from Combo_Base import Combo
from DynamicModel_Package.DynamicModel2D_Base import DynamicModel2D
import sympy as sp

x,y,r,theta,dr,dtheta = sp.symbols('x y r theta dr dtheta')

class Combo_plot ():
    def __init__ (self,combo:Combo,parameters_dict:dict = {}):
        self.combo = combo
        self.parameters_dict = parameters_dict
        self.model = self.create_model()
        self.parameter_values = list(self.parameters_dict.values())
        self.convert_nullclines()

    def deriv_x (self, vars:dict,pars:dict):
        if len(pars)>0:
            dict_all = vars.copy()
            dict_all.update(pars)
        else:
            dict_all = vars
        expr_str = str(self.combo.P.as_expr())
        expr = sp.sympify(expr_str)
        ret = expr.subs(dict_all)
        return ret
    
    def deriv_y (self, vars:dict,pars:dict):
        if len(pars)>0:
            dict_all = vars.copy()
            dict_all.update(pars)
        else:
            dict_all = vars
        expr_str = str(self.combo.Q.as_expr())
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
    


