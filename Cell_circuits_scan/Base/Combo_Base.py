import numpy as np
import sympy as sp
from coeff_base import coeff

x,y,r,theta,dr,dtheta = sp.symbols('x y r theta dr dtheta',real = True)

class Combo ():
    _id_counter = 1  # Class variable to keep track of the next ID

    def __init__ (self,ps,id_counter=None):
        self.ps = ps
        self.num_parameters = len(self.ps)
        self.P, self.Q, self.variables, self.parameters = self.recreate_polynomial()
        self.polys = [self.P,self.Q]
        self.ps_list = [p.p for p in self.ps]
        self.ps_list_no_symbol = [p.rest for p in ps]
        self.ps_symbols_list = [p.symbol for p in ps]
        self.ps_letters = [p.letter for p in ps]
        self.dr = dr
        self.dtheta = dtheta
        self.rdot = None
        self.thetadot = None
        self.jac = None
        self.nullclines_x = None
        self.nullclines_y = None
        self.x = x
        self.y = y
        self.r = r
        self.theta = theta
        self.id = self.create_id()
        if id_counter==None:
            self.id_counter = Combo._id_counter  # Assign the current value of the counter to the instance
        else:
            self.id_counter = id_counter
        Combo._id_counter += 1
    

    @classmethod
    def reset_id_counter(cls):
        cls._id_counter = 1  # Reset the counter to 1

    def to_list_strings(self):
        return [p.p for p in self.ps]
    
    def create_id (self):
        #creates a string that is unique for each combo
        #each place return the sign of the corresponding element +1
        #the elements are in the following order: [a00,a10,a20,a01,a02,a11,b00,b10,b20,b01,b02,b11]
        elements = ['a00','a10','a20','a01','a02','a11','b00','b10','b20','b01','b02','b11']
        id = ""
        for element in elements:
            if element in self.ps_list_no_symbol:
                symbol = self.ps_symbols_list[self.ps_list_no_symbol.index(element)]
                if symbol=='-':
                    id+=str(0)
                else:
                    id+=str(2)
            else:
                id+=str(1)
        return id

    def recreate_polynomial (self):
        P = 0 # dxdt polynomial
        Q = 0 # dydt polynomial
        # Define the polynomials from the paramters:
        variables = [x,y]
        params = []
        for i,p in enumerate(self.ps):
            num1 = int(p.number1)
            num2 = int(p.number2)
            if i<3:
                if p.letter=='a':
                    if p.symbol == '+':
                        P += (x**(num1+1))*(y**num2)
                    else:
                        P -= (x**(num1+1))*(y**num2)
                else:
                    if p.symbol == '+':
                        Q += (x**num1)*(y**(num2+1))
                    else:
                        Q -= (x**num1)*(y**(num2+1))
            else:
                dp = sp.symbols(f'dp{i-3}',positive = True)
                params.append(dp)
                if p.letter=='a':
                    if p.symbol == '+':
                        P += (dp*x**(num1+1))*(y**num2)
                    else:
                        P -= (dp*x**(num1+1))*(y**num2)
                else:
                    if p.symbol == '+':
                        Q += (dp*x**num1)*(y**(num2+1))
                    else:
                        Q -= (dp*x**num1)*(y**(num2+1))
        P = sp.Poly(P,x,y)
        Q = sp.Poly(Q,x,y)    
        return P,Q,variables,params

    def poly_cart_to_polar (self,P = None,Q = None):
        flag_return = True
        if P == None:
            P = self.P
            flag_return = False
        if Q == None:
            Q = self.Q
            flag_return=False
        P = P.as_expr()
        Q = Q.as_expr()

        # Define the transformations
        x_polar = r * sp.cos(theta)
        y_polar = r * sp.sin(theta)

        # Substitute x and y with their polar forms in P3 and Q3
        P_polar = P.subs({x: x_polar, y: y_polar})
        Q_polar = Q.subs({x: x_polar, y: y_polar})

        # Compute the derivatives of the transformations
        dx_polar = sp.diff(x_polar, r) * dr + sp.diff(x_polar, theta) * dtheta
        dy_polar = sp.diff(y_polar, r) * dr + sp.diff(y_polar, theta) * dtheta

        eq1 = sp.Eq(P_polar,dx_polar)
        eq2 = sp.Eq(Q_polar,dy_polar)
        solution = sp.solve((eq1, eq2), (dr, dtheta))
        if flag_return == False:
            self.rdot = solution[dr]
            self.thetadot = solution[dtheta]
        return solution

    def highest_power_polar (self,P = None,Q = None):
        flag_use_self = False
        if P==None:
            P = self.P
            flag_use_self = True
        if Q==None:
            Q = self.Q
        if flag_use_self:
            if self.rdot == None:
                self.poly_cart_to_polar()
            rdot = self.rdot
        else:
            sol = self.poly_cart_to_polar(P,Q)
            rdot = sol[dr]
        coeffs = sp.Poly(rdot,r).all_coeffs()
        highest_power = len(coeffs)-1
        highest_power_coeff = coeffs[0]
        highest_power_coeff = sp.trigsimp(highest_power_coeff)
        return highest_power_coeff,highest_power
    
    def first_nonzero_index(self,lst):
        for index, value in enumerate(lst):
            if value != 0:
                return index
        return None  # Return None if all elements are zero

    def lowest_power_polar (self,P = None,Q = None):
        flag_use_self = False
        if P==None:
            P = self.P
            flag_use_self = True
        if Q==None:
            Q = self.Q
        if flag_use_self:
            if self.rdot == None:
                self.poly_cart_to_polar()
            rdot = self.rdot
        else:
            sol = self.poly_cart_to_polar(P,Q)
            rdot = sol[dr]
        coeffs = sp.Poly(rdot,r).all_coeffs()
        coeffs_ascending = coeffs[::-1]
        lowest_power = self.first_nonzero_index(coeffs_ascending)
        lowest_power_coeff = coeffs_ascending[lowest_power]
        lowest_power_coeff = sp.trigsimp(lowest_power_coeff)
        return lowest_power_coeff,lowest_power

    def create_jac (self):
        P = self.P.as_expr()
        Q = self.Q.as_expr()
        polys_ = sp.Matrix([P,Q])
        jac = polys_.jacobian(self.variables)
        self.jac = jac

    def check_fixed_point_stable(self,fp):
        if self.jac==None:
            self.create_jac()
        jac_subs = self.jac.subs({x:fp[0].as_expr(),y:fp[1].as_expr()})
        evs = jac_subs.eigenvals(jac_subs)
        print(evs)
        ret_stable = True
        for ev in evs:
            if ev.is_positive:
                ret_stable = False
        return ret_stable

    def find_nullclines (self):
        nc_x_x = sp.solve(sp.Eq(self.P,0),x) #result as a function of y
        nc_x_y = sp.solve(sp.Eq(self.P,0),y) #result as a function of x
        nc_y_x = sp.solve(sp.Eq(self.Q,0),x)
        nc_y_y = sp.solve(sp.Eq(self.Q,0),y)
        nc_x = []
        nc_y = []
        for n in nc_x_x:
            nc_x.append((n,y))
        for n in nc_x_y:
            nc_x.append((x,n))
        for n in nc_y_x:
            nc_y.append((n,y))
        for n in nc_y_y:
            nc_y.append((x,n))
        self.nullclines_x = nc_x
        self.nullclines_y = nc_y
        nullclines = []
        for n in nc_x: nullclines.append(n)
        for n in nc_y: nullclines.append(n)
        self.nullclines = nullclines

