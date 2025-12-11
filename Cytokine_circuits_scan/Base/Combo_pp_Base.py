import numpy as np
import sympy as sp
import time
from multiprocessing import Pool, TimeoutError as MPTimeoutError

x,y = sp.symbols('x y')

class Combo_pp ():
    _id_counter = 1  # Class variable to keep track of the next ID

    def __init__ (self,xdot,ydot,
                  parameters,
                  num_interactions,
                  x = sp.symbols('x'),
                  y = sp.symbols('y'),
                  parameter_placement = 'def',
                  find_nullclines = False):
        self.x = x
        self.y = y
        self.variables = [x,y]
        self.xdot = xdot
        self.ydot = ydot
        self.parameters = parameters
        self.num_dof = num_interactions
        self.parameter_placement = parameter_placement
        if find_nullclines:
            self.find_nullclines()
        self.id_counter = Combo_pp._id_counter
        Combo_pp._id_counter += 1

    @classmethod
    def reset_id_counter(cls):
        cls._id_counter = 1  # Reset the counter to 1

    def find_nullclines (self):
        #nc_x_x = sp.solve(sp.Eq(self.xdot,0),x) #result as a function of y
        nc_x_y = sp.solve(sp.Eq(self.xdot,0),y) #result as a function of x
        #nc_y_x = sp.solve(sp.Eq(self.ydot,0),x)
        nc_y_y = sp.solve(sp.Eq(self.ydot,0),y)
        nc_x = []
        nc_y = []
        #for n in nc_x_x:
            #nc_x.append((n,y))
        for n in nc_x_y:
            nc_x.append((x,n))
        #for n in nc_y_x:
            #nc_y.append((n,y))
        for n in nc_y_y:
            nc_y.append((x,n))
        self.nullclines_x = nc_x
        self.nullclines_y = nc_y
        nullclines = []
        for n in nc_x: nullclines.append(n)
        for n in nc_y: nullclines.append(n)
        self.nullclines = nullclines
