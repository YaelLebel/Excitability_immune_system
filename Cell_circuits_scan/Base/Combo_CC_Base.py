import numpy as np
import sympy as sp
import time
from multiprocessing import Pool, TimeoutError as MPTimeoutError

x,y = sp.symbols('x y')

class Combo_CC ():
    _id_counter = 1  # Class variable to keep track of the next ID

    def __init__ (self,prod_x_terms,rem_x_terms,
                  prod_y_terms,rem_y_terms,
                  CC_x=100,CC_y=100,
                  parameter_placement = 'y',
                  find_nullclines = True):
        self.x = x
        self.y = y
        self.variables = [x,y]
        self.prod_x_terms = prod_x_terms
        self.rem_x_terms = rem_x_terms
        self.prod_y_terms = prod_y_terms
        self.rem_y_terms = rem_y_terms
        self.C_x = CC_x
        self.C_y = CC_y
        self.num_dof = len(prod_x_terms)+len(prod_y_terms)+len(rem_x_terms)+len(rem_y_terms)
        self.parameter_placement = parameter_placement
        self.P,self.Q,self.params = self.create_polynomials()
        if find_nullclines:
            self.find_nullclines()
        self.id_counter = Combo_CC._id_counter
        Combo_CC._id_counter += 1

    @classmethod
    def reset_id_counter(cls):
        cls._id_counter = 1  # Reset the counter to 1

    def create_polynomials(self):
        if self.num_dof<=3:
            return self.create_polynomials_no_parameters()
        elif self.parameter_placement=='y':
            return self.create_polynomials_parameter_placement_y()
        else:
            return self.create_polynomials_parameter_placement_x()

    def create_polynomials_no_parameters (self):
        params = []
        P = 0
        Q = 0
        for term in self.prod_x_terms:
            term = sp.sympify(term)
            if self.C_x==None:
                P += term
            else:
                P += term*(1-x/self.C_x)
        for term in self.rem_x_terms:
            term = sp.sympify(term)
            P -= term
        P = x*P
        for term in self.prod_y_terms:
            term = sp.sympify(term)
            if self.C_y==None:
                Q += term
            else:
                Q += term*(1-y/self.C_y)
        for term in self.rem_y_terms:
            term = sp.sympify(term)
            Q -= term
        Q = y*Q
        return P,Q,params

    def create_polynomials_parameter_placement_x (self):
        pars = sp.symbols(('dp0:%d'%(self.num_dof-3)),positive = True)
        if self.num_dof==1:
            pars = list(pars)
        P = 0
        Q = 0
        counter=0
        prod_x = sp.sympify(0)
        for term in self.prod_x_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                prod_x += dp*term
            else:
                prod_x += term
            counter += 1
        rem_x = sp.sympify(0)
        for term in self.rem_x_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                rem_x += dp*term
            else:
                rem_x += term
            counter += 1
        if self.C_x != None:
            P += prod_x*(1-x/self.C_x)
        else:
            P += prod_x
        P -= rem_x
        P = x*P
        prod_y = sp.sympify(0)
        for term in self.prod_y_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                prod_y += dp*term
            else:
                prod_y += term
            counter += 1
        rem_y = sp.sympify(0)
        for term in self.rem_y_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                rem_y += dp*term
            else:
                rem_y += term
            counter += 1
        if self.C_y != None:
            Q += prod_y*(1-y/self.C_y)
        else:
            Q += prod_y
        Q -= rem_y
        Q = y*Q
        return P,Q,pars

    def create_polynomials_parameter_placement_y (self):
        pars = sp.symbols(('dp0:%d'%(self.num_dof-3)),positive = True)
        if self.num_dof==1:
            pars = list(pars)
        P = 0
        Q = 0
        counter=0
        prod_y = sp.sympify(0)
        for term in self.prod_y_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                prod_y += dp*term
            else:
                prod_y += term
            counter += 1
        rem_y = sp.sympify(0)
        for term in self.rem_y_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                rem_y += dp*term
            else:
                rem_y += term
            counter += 1
        if self.C_y != None:
            Q += prod_y*(1-y/self.C_y)
        else:
            Q += prod_y
        Q -= rem_y
        Q = y*Q
        prod_x = sp.sympify(0)
        for term in self.prod_x_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                prod_x += dp*term
            else:
                prod_x += term
            counter += 1
        rem_x = sp.sympify(0)
        for term in self.rem_x_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                rem_x += dp*term
            else:
                rem_x += term
            counter += 1
        if self.C_x != None:
            P += prod_x*(1-x/self.C_x)
        else:
            P += prod_x
        P -= rem_x
        P = x*P
        return P,Q,pars

    def create_polynomials_near_origin (self):
        if self.num_dof<=3:
            return self.create_polynomials_near_origin_no_parameters()
        elif self.parameter_placement=='y':
            return self.create_polynomials_near_origin_parameter_placement_y()
        else:
            return self.create_polynomials_near_origin_parameter_placement_x()

    def create_polynomials_near_origin_no_parameters (self):
        params = []
        P = 0
        Q = 0
        for term in self.prod_x_terms:
            term = sp.sympify(term)
            P += term
        for term in self.rem_x_terms:
            term = sp.sympify(term)
            P -= term
        P = x*P
        for term in self.prod_y_terms:
            term = sp.sympify(term)
            Q += term
        for term in self.rem_y_terms:
            term = sp.sympify(term)
            Q -= term
        Q = y*Q
        return P,Q,params

    def create_polynomials_near_origin_parameter_placement_x (self):
        pars = sp.symbols(('dp0:%d'%(self.num_dof-3)),positive = True)
        if self.num_dof==1:
            pars = list(pars)
        P = 0
        Q = 0
        counter=0
        prod_x = sp.sympify(0)
        for term in self.prod_x_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                prod_x += dp*term
            else:
                prod_x += term
            counter += 1
        rem_x = sp.sympify(0)
        for term in self.rem_x_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                rem_x += dp*term
            else:
                rem_x += term
            counter += 1
        P += prod_x
        P -= rem_x
        P = x*P
        prod_y = sp.sympify(0)
        for term in self.prod_y_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                prod_y += dp*term
            else:
                prod_y += term
            counter += 1
        rem_y = sp.sympify(0)
        for term in self.rem_y_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                rem_y += dp*term
            else:
                rem_y += term
            counter += 1
        Q += prod_y
        Q -= rem_y
        Q = y*Q
        return P,Q,pars

    def create_polynomials_near_origin_parameter_placement_y (self):
        pars = sp.symbols(('dp0:%d'%(self.num_dof-3)),positive = True)
        if self.num_dof==1:
            pars = list(pars)
        P = 0
        Q = 0
        counter=0
        prod_y = sp.sympify(0)
        for term in self.prod_y_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                prod_y += dp*term
            else:
                prod_y += term
            counter += 1
        rem_y = sp.sympify(0)
        for term in self.rem_y_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                rem_y += dp*term
            else:
                rem_y += term
            counter += 1
        Q += prod_y
        Q -= rem_y
        Q = y*Q
        prod_x = sp.sympify(0)
        for term in self.prod_x_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                prod_x += dp*term
            else:
                prod_x += term
            counter += 1
        rem_x = sp.sympify(0)
        for term in self.rem_x_terms:
            term = sp.sympify(term)
            if counter<len(pars):
                dp = pars[counter]
                rem_x += dp*term
            else:
                rem_x += term
            counter += 1
        P += prod_x
        P -= rem_x
        P = x*P
        return P,Q,pars

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

    def create_jac (self):
        P = self.P.as_expr()
        Q = self.Q.as_expr()
        polys_ = sp.Matrix([P,Q])
        jac = polys_.jacobian(self.variables)
        self.jac = jac

    def linear_stability (self,fp,f=None,g=None):
        if f==None:
            f = self.P.as_expr()
        if g==None:
            g = self.Q.as_expr()
        J = sp.Matrix([[sp.diff(f, x), sp.diff(f, y)],
                    [sp.diff(g, x), sp.diff(g, y)]])
        J_origin = J.subs({x:fp[0] if (type(fp[0])==int or type(fp[0])==float) else fp[0].as_expr(),y:fp[1] if (type(fp[1])==int or type(fp[1])==float) else fp[1].as_expr()})
        eigenvals = J_origin.eigenvals()
        eigenvals_reals = [sp.re(x) for x in eigenvals.keys()]
        if 0 in eigenvals_reals:
            return None
        else:
            if np.all([x<0 for x in eigenvals_reals]):
                return True
            else:
                return False

    def center_manifold_stability(self,fp,f=None,g=None,max_level = 10):
        if f==None:
            f = self.P.as_expr()
        if g==None:
            g = self.Q.as_expr()

        J = sp.Matrix([[sp.diff(f, x), sp.diff(f, y)],
                    [sp.diff(g, x), sp.diff(g, y)]])

        J_origin = J.subs({x:fp[0] if (type(fp[0])==int or type(fp[0])==float) else fp[0].as_expr(),y:fp[1] if (type(fp[1])==int or type(fp[1])==float) else fp[1].as_expr()})

        eigenvals = J_origin.eigenvals()

        A = list(eigenvals.keys())[0]
        B = list(eigenvals.keys())[1] if len(eigenvals)>1 else list(eigenvals.keys())[0]
        flag = True
        level = 2
        time_limit = 20

        while flag and level<=max_level:
            print("checking order ",level)
            pars = sp.symbols(('h2:%d'%(level+1)))
            if level<3:
                pars = list(pars)
            h_x = 0
            for i,par in enumerate(pars):
                h_x += par*x**(i+2)
            h_prime_x = sp.diff(h_x,x)
            f_reduced = sp.Poly(f.subs(y,h_x),x)
            g_reduced = sp.Poly(g.subs(y,h_x),x)
            Nhx = h_prime_x*(A*x+f_reduced)-(B*h_x+g_reduced)
            coeffs = Nhx.all_coeffs()
            Nhx_reduced = sp.sympify(0)
            i=0
            while pars[-1] not in Nhx_reduced.free_symbols:
                Nhx_reduced += coeffs[::-1][i]*x**(i+2)
                i+=1

            Nhx_reduced = sp.Poly(Nhx_reduced,x)
            coeffs = Nhx_reduced.coeffs()
            sol = sp.solve(coeffs,pars)
            h_solvd = 0
            if type(sol)==dict:
                for i,par in enumerate(pars):
                    h_solvd += sol[par]*x**(i+2)
            else:
                for i in range(len(sol[0])):
                    h_solvd += sol[0][i]*x**(i+2)
            if h_solvd==0:
                level +=1
            else:
                flag = False
        f_center = sp.Poly(f.subs(y,h_solvd),x)
        print("center manifold approximation ",f_center)
        stability = sp.sign(f_center.subs(x,0.001))
        if stability==0:
            return None
        else: return not stability==1

    def check_stability(self,fp,max_level = 10):
        if fp[0]==0 and fp[1]==0:
            f,g,_ = self.create_polynomials_near_origin()
        else:
            f = self.P.as_expr()
            g = self.Q.as_expr()
        lin_stab = self.linear_stability(fp,f=f,g=g)
        if lin_stab==None:
            cen_man_stab = self.center_manifold_stability(fp,f=f,g=g,max_level=max_level)
            return cen_man_stab
        else:
            return lin_stab
