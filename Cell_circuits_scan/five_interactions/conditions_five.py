import sympy as sp
import sys
import os
import time
import timeout_decorator
import numpy as np
from itertools import product

# Get the absolute path to the 'Base' directory
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Base'))

# Add the 'Base' directory to sys.path
sys.path.append(base_dir)

def linear_stability (cm,fp,f,g,verbose = False):
    if f==None:
        f = cm.P.as_expr()
    if g==None:
        g = cm.Q.as_expr()
    x = cm.x
    y = cm.y
    J = sp.Matrix([[sp.diff(f, x), sp.diff(f, y)],
                [sp.diff(g, x), sp.diff(g, y)]])
    J_origin = J.subs({x:fp[0] if (type(fp[0])==int or type(fp[0])==float) else fp[0].as_expr(),y:fp[1] if (type(fp[1])==int or type(fp[1])==float) else fp[1].as_expr()})
    eigenvals = J_origin.diagonal()
    eigenvals_reals = [sp.re(x) for x in eigenvals]
    print("eigenvalues real part:",eigenvals_reals)
    if verbose:
        print("jacobian at origin:",J_origin)
        print("eigenvalues:",eigenvals)
    if np.any([sp.sign(ev)==1 for ev in eigenvals_reals]):
        return False
    if np.any([ev==0 for ev in eigenvals_reals]):
        return None
    if np.all([sp.sign(ev)==-1 for ev in eigenvals_reals]):
        return True
    return eigenvals_reals



def check_stability_only_decreasing_y (cm,fp,verbose = False,
                                      check_x_greater_than_point = True,
                                      epsilon = 0.01,check_linear_approximation = True):
    f = cm.P.as_expr()
    g = cm.Q.as_expr()
    if check_linear_approximation:  
        lin_stab = linear_stability(cm,[0,0],f=f,g=g,verbose=verbose)
    if (check_linear_approximation==False) or (lin_stab==None):
        print("linear stability failed or skipped. Checking x behaviour explicitly")
        if verbose:
            print ("xdot near origin:",sp.diff(f,cm.x).subs(cm.y,fp[1]))
        xdot_near_fp = sp.diff(f,cm.x).subs(cm.y,fp[1]).subs(cm.x,fp[0]+(2*int(check_x_greater_than_point)-1)*epsilon)
        if verbose:
            print("estimating xdot at:",fp[0]+(2*int(check_x_greater_than_point)-1)*epsilon)
            print("estimated xdot at point proximity:",xdot_near_fp)
        if len(xdot_near_fp.free_symbols)==0:
            if sp.sign(xdot_near_fp)==0:
                return 0
            else:
                if check_x_greater_than_point:
                    return sp.sign(xdot_near_fp)==-1
                else:
                    return sp.sign(xdot_near_fp)==1
        else:
            print("stability dependant of parameter values; checking different parameter values")
            if cm.C_x!=None:
                par_vals = [0.001,1.01,0.99*cm.C_x,100*cm.C_x]
            elif cm.C_y!=None:
                par_vals = [0.001,1.01,0.99*cm.C_y,100*cm.C_y]
            else:
                par_vals = [0.001,1.01,0.99*100,100*100]
            combinations = product(par_vals, repeat=len(xdot_near_fp.free_symbols))
            dict_list = [dict(zip(xdot_near_fp.free_symbols, combo)) for combo in combinations]
            xdot_near_fp_vals = [xdot_near_fp.subs(d).subs(cm.x,fp[0]+(2*int(check_x_greater_than_point)-1)*epsilon) for d in dict_list]
            print("xdot near origin on different parameter values:",xdot_near_fp_vals)
            if len(set([sp.sign(f_cen_val) for f_cen_val in xdot_near_fp_vals]))==1:
                if check_x_greater_than_point:
                    return sp.sign(xdot_near_fp)==-1
                else:
                    return sp.sign(xdot_near_fp)==1
            else:
                return None
    else:
        print("stability determined by linear approximation:",lin_stab)
        return lin_stab

def center_manifold_stability_level (cm,level,A,B,verbose = False):
    x = cm.x
    y = cm.y
    xdot = cm.P.as_expr()
    ydot = cm.Q.as_expr()
    pars = sp.symbols(('h2:%d'%(level+1)))
    if level<3:
        pars = list(pars)
    h = 0
    if B==0:
        #construct y = h(x)
        for i,par in enumerate(pars):
            h += par*x**(i+2)
        h_tag = sp.diff(h,x)
        Nhx = ydot.subs(y,h_tag)-h_tag*xdot.subs(y,h_tag)
        Nhx = sp.Poly(Nhx,x)
        coeffs = Nhx.coeffs()
        if verbose:
            print("all coeffs:",coeffs)
        i=2
        h_solvd = 0
        while i<len(coeffs) and h_solvd==0:
            solutions = sp.solve(Nhx.coeffs()[::-1][:i],pars[:i])
            if verbose:
                print("coeffs subset:",Nhx.coeffs()[::-1][:i])
            if verbose:
                print("solution for this subset:",solutions)
            if type(solutions)==list:
                if type(solutions[0])==dict:
                    for d in solutions:
                        if np.any([v!=0 for v in d.values()]):
                            #check if this solution is non-zero (at least one of the values)
                            for power,v in enumerate(d.values()):
                                h_solvd += v*x**(power+2)
                            break
                elif type(solutions[0])==tuple:
                    for s in solutions:
                        if np.any([v!=0 for v in s]):
                            for power,v in enumerate(s):
                                h_solvd += v*x**(power+2)
                            break
            elif type(solutions)==dict:
                if np.any([v!=0 for v in solutions.values()]):
                    for power,v in enumerate(solutions.values()):
                        h_solvd += v*x**(power+2)
            elif type(solutions)==tuple:
                if np.any([v!=0 for v in solutions]):
                    for power,v in enumerate(s):
                        h_solvd += v*x**(power+2)
            if verbose:
                print("h:",h_solvd)

            i+=1

    if A==0:
        #construct x = h(y)
        for i,par in enumerate(pars):
            h += par*y**(i+2)
        h_tag = sp.diff(h,y)
        Nhy = xdot.subs(x,h_tag)-h_tag*ydot.subs(x,h_tag)
        Nhy = sp.Poly(Nhy,y)
        coeffs = Nhy.coeffs()
        i=1
        h_solvd = 0
        if verbose:
            print("all coeffs:",coeffs)
        while i<len(coeffs) and h_solvd==0:
            solutions = sp.solve(Nhy.coeffs()[::-1][:i],pars[:i])
            if verbose:
                print("coeffs subset:",Nhy.coeffs()[::-1][:i])
            if verbose:
                print("solution for this subset:",solutions)
            if type(solutions)==list:
                if type(solutions[0])==dict:
                    for d in solutions:
                        if np.any([v!=0 for v in d.values()]):
                            for power,v in enumerate(d.values()):
                                h_solvd += v*y**(power+2)
                elif type(solutions[0])==tuple:
                    for s in solutions:
                        if np.any([v!=0 for v in s]):
                            for power,v in enumerate(s):
                                h_solvd += v*y**(power+2)
            elif type(solutions)==dict:
                if np.any([v!=0 for v in solutions.values()]):
                    for power,v in enumerate(solutions.values()):
                        h_solvd += v*y**(power+2)
            elif type(solutions)==tuple:
                if np.any([v!=0 for v in solutions]):
                    for power,v in enumerate(s):
                        h_solvd += v*y**(power+2)
            i+=1

    return h_solvd

def get_lowest_approx_poly (pol,var):
    lowest_order = next((i for i, x in enumerate(pol.all_coeffs()[::-1]) if x), None)
    pol_lowest_coeff= pol.all_coeffs()[::-1][lowest_order]*var**lowest_order
    return pol_lowest_coeff

@timeout_decorator.timeout(5)
def center_manifold_stability(cm,fp,f=None,g=None,max_level = 10,verbose = False):
    if f==None:
        f = cm.P.as_expr()
    if g==None:
        g = cm.Q.as_expr()
    x = cm.x
    y = cm.y

    J = sp.Matrix([[sp.diff(f, x), sp.diff(f, y)],
                [sp.diff(g, x), sp.diff(g, y)]])

    J_origin = J.subs({x:fp[0] if (type(fp[0])==int or type(fp[0])==float) else fp[0].as_expr(),y:fp[1] if (type(fp[1])==int or type(fp[1])==float) else fp[1].as_expr()})

    eigenvals = J_origin.diagonal()

    A = eigenvals[0]
    B = eigenvals[1] if len(eigenvals)>1 else eigenvals[0]
    if verbose:
        print(f"A:{A},B:{B}")
    flag = True
    level = 2
    time_limit = 3
    while flag and level<=max_level:
        try:
            h_solvd = center_manifold_stability_level(cm,level,A,B,verbose)
            if verbose:
                print(f"h_solvd main function, level {level}: {h_solvd}")
            if h_solvd ==0:
                level +=1
            else:
                flag = False
                level +=1
        except TimeoutError:
            h_solvd=0
            flag = False
    if A==0:
        f_center = sp.Poly(f.subs(x,h_solvd),y)
        if verbose:
            print("center manifold approximation:",h_solvd)
            print("center manifold reduction ",f_center)

        if len(f_center.subs(y,0.001).free_symbols)==0:
            if sp.sign(f_center.subs(y,0.001))==0:
                return 0
            else:
                return sp.sign(f_center.subs(y,0.001))==-1
        else:
            f_center = get_lowest_approx_poly(f_center,y)
            #will check the sign for different values of dp0: very small, order of 1, order of carrying capacity, larger than carrying capacity
            print("stability dependant of dp0 values; checking different dp0 values")
            if cm.C_y!=None:
                par_vals = [0.001,1.01,0.99*cm.C_y,100*cm.C_y]
            elif cm.C_x!=None:
                par_vals = [0.001,1.01,0.99*cm.C_x,100*cm.C_x]
            else:
                par_vals = [0.001,1.01,0.99*100,100*100]
            combinations = product(par_vals, repeat=len(f_center.free_symbols))
            dict_list = [dict(zip(f_center.free_symbols, combo)) for combo in combinations]
            f_center_vals = [f_center.subs(d).subs(cm.x,fp[0]+0.001) for d in dict_list]
            if verbose:
                print("center manifold on different parameter values:",f_center_vals)
            if len(set([sp.sign(f_cen_val) for f_cen_val in f_center_vals]))==1:
                return sp.sign(f_center_vals[0])==-1
            else:
                return None
    if B==0:
        g_center = sp.Poly(g.subs(y,h_solvd),x)
        if verbose:
            print("center manifold approximation:",h_solvd)
            print("center manifold reduction ",g_center)
        if len(g_center.subs(x,0.001).free_symbols)==0:
            if sp.sign(g_center.subs(x,0.001))==0:
                return 0
            else:
                return sp.sign(g_center.subs(x,0.001))==-1
        else:
            g_center = get_lowest_approx_poly(g_center,x)
            #will check the sign for different values of dp0: very small, order of 1, order of carrying capacity, larger than carrying capacity
            print("stability dependant of dp0 values; checking different dp0 values")
            if cm.C_x!=None:
                par_vals = [0.001,1.01,0.99*cm.C_x,100*cm.C_x]
            elif cm.C_y!=None:
                par_vals = [0.001,1.01,0.99*cm.C_y,100*cm.C_y]
            else:
                par_vals = [0.001,1.01,0.99*100,100*100]
            combinations = product(par_vals, repeat=len(g_center.free_symbols))
            dict_list = [dict(zip(g_center.free_symbols, combo)) for combo in combinations]
            g_center_vals = [g_center.subs(d).subs(cm.y,fp[1]+0.001) for d in dict_list]
            if verbose:
                print("center manifold on different parameter values:",g_center_vals)
            if len(set([sp.sign(g_cen_val) for g_cen_val in g_center_vals]))==1:
                return sp.sign(g_center_vals[0])==-1
            else:
                return None

def check_stability(cm,fp,max_level = 10,verbose = False):
    if fp[0]==0 and fp[1]==0:
        f,g,_ = cm.create_polynomials_near_origin()
    else:
        f = cm.P.as_expr()
        g = cm.Q.as_expr()
    lin_stab = linear_stability(cm,fp,f=f,g=g,verbose=verbose)
    if lin_stab==None:
        print("linear stability failed. Checking center manifold")
        try:
            cen_man_stab = center_manifold_stability(cm,fp,f=f,g=g,max_level=max_level,verbose=verbose)
            return cen_man_stab
        except Exception as e:
            print("exception:",e)
            return None
    else:
        print("stability determined by linear approximation:",lin_stab)
        return lin_stab

