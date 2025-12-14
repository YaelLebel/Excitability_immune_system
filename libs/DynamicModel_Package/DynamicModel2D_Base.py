from DynamicModel_Package.ModelBase import DynamicModel
import matplotlib.pyplot as plt
import numpy as np

class DynamicModel2D (DynamicModel):
    def __init__ (self,variable_x:str,derivative_function_x:"function",parameters_x:dict,
                  variable_y:str,derivative_function_y:"function",parameters_y:dict):
        DynamicModel.__init__(self)
        self.add_variable(variable_x,derivative_function_x,parameters_x,'x')
        self.add_variable(variable_y,derivative_function_y,parameters_y,'y')
        self.variable_x = variable_x
        self.variable_y = variable_y
        self.nullclines = {}
        self.parameters_x = parameters_x
        self.parameters_y = parameters_y
    
    def add_variable(self, name:str, derivative_function:"function", parameters:dict,axis:str):
        """
        Add a variable to the model with its time derivative function and parameters.
        
        Parameters:
            name (str): Name of the variable.
            derivative_function (function): Function describing the time derivative of the variable.
            parameters (dict): Parameters required by the derivative function.
            axis(str): axis that belongs to this variable (x or y)
        """
        self.variables[name] = {
            'derivative_function': derivative_function,
            'parameters': parameters,
            'axis':axis
        }
    
    def add_nullcline(self, name:str, nullcline_function:"function", parameters:dict,y_fun_x:bool):
        """
        Add a nullcline of a variable to the model with its function and parameters.
        
        Parameters:
            name (str): Name of the variable.
            nullcline_function (function): Function describing the nullcline of the variable.
            parameters (dict): Parameters required by the derivative function.
            y_fun_x (Boolean) : wether the function is provided as y(x), according to the axis defined for each variable
        """
        
        mydict = {
            'nullcline_function': nullcline_function,
            'parameters': parameters,
            'y_fun_x': y_fun_x
        }    
        if name in self.nullclines.keys():
            self.nullclines[name].append(mydict)
        else:
            self.nullclines[name] = [mydict]
    
    def create_meshgrid_derivatives (self,n_X:int,n_Y:int,Xlim:tuple,Ylim:tuple,t:float):
        """
        Parameters
        ----------
        n_X : int
            number of points to evaluate in p axis.
        n_Y : int
            number of points to evaluate in  H axis.
        Xlim : tuple (float,float)
            limits to evaluate X.
        Ylim : tuple (float,float)
            limits to evaluate Y.


        Returns
        -------
        Xv, Yv : np.array (np.meshgrid)
        Xdot,Ydot : np.array (np.meshgrid - like)
            derivatives in each point 
        """
        Xv,Hv = np.meshgrid(np.linspace(Xlim[0],Xlim[1],n_X),np.linspace(Ylim[0],Ylim[1],n_Y))
        Xdot = np.empty_like(Xv)
        Ydot = np.empty_like(Hv)

        X_name = self.variable_x
        Y_name = self.variable_y
        for i in range(n_X):
            for j in range(n_Y):
                var_dict = {'t':t,X_name:Xv[i,j],Y_name:Hv[i,j]}
                derivatives = super().compute_derivative(var_dict)
                Ydot[i,j] = derivatives[Y_name]
                Xdot[i,j] = derivatives[X_name]
        return Xv,Hv,Xdot,Ydot
    
    
    def plot_nullclines (self,variable:str,xlim:tuple,ylim:tuple, fig = None, ax= None, color = None,plot_args = None):
        flag_ret_figs = False
        if fig == None and ax == None:
            fig,ax = plt.subplots()
            flag_ret_figs = True
        
        for i,nullcline_dict in enumerate(self.nullclines[variable]):
            f = nullcline_dict['nullcline_function']
            pars = nullcline_dict['parameters']

            y_function_x = nullcline_dict['y_fun_x']
            if y_function_x:
                x = np.linspace(min(xlim[0],ylim[0]),max(xlim[1],ylim[1]),1000)
                y = [f(_,pars) for _ in x]
                mask = np.isreal(y)
                x = x[mask]
                y = np.array(y)[mask]
                import sympy as sp
                x_to_keep = []
                y_to_keep = []
                for x_i,y_i in zip(x,y):
                    if sp.im(y_i)==0:
                        x_to_keep.append(x_i)
                        y_to_keep.append(y_i)
                x = np.array(x_to_keep)
                y = np.array(y_to_keep)
                if color != None:
                    if i==0:
                        ax.plot(x,y,color = color,label = f'{variable} nullcline',zorder = 1,**plot_args) 
                    else:
                        ax.plot(x,y,color = color,zorder = 1,**plot_args)
                else:
                    ax.plot(x,y,label = f'{variable} nullcline',zorder = 1,**plot_args)
            else:
                x = np.linspace(min(xlim[0],ylim[0]),max(xlim[1],ylim[1]),1000)
                y = [f(_,pars) for _ in x]
                mask = np.isreal(y)
                x = x[mask]
                y = np.array(y)[mask]

                import sympy as sp
                x_to_keep = []
                y_to_keep = []
                for x_i,y_i in zip(x,y):
                    if sp.im(y_i)==0:
                        x_to_keep.append(x_i)
                        y_to_keep.append(y_i)
                x = np.array(x_to_keep)
                y = np.array(y_to_keep)
                if color != None:
                    if i==0:
                        ax.plot(y,x, color = color,label = f'{variable} nullcline',zorder = 1,**plot_args)
                    else:
                        ax.plot(y,x, color = color,zorder = 1,**plot_args)
                else:
                    ax.plot(y,x,label = f'{variable} nullcline',zorder = 1,**plot_args)
        ax.legend(loc = 'upper right')
        if flag_ret_figs:
            return fig, ax
        
    
    def plot_streamplot (self, t:float,n_X:int,n_Y:int, xlim:tuple, ylim:tuple, fig=None, ax=None,density = 1,alpha = 1,color = 'grey'):
        xv,yv,xdot,ydot = self.create_meshgrid_derivatives(n_X,n_Y,xlim,ylim,t)        
        flag_ret_figs = False
        if fig == None and ax == None:
            fig,ax = plt.subplots()
            flag_ret_figs = True
        
        ax.streamplot(xv,yv,xdot,ydot,color = color,density = density,zorder = 0,)
        
        if flag_ret_figs:
            return fig, ax
    
    
    def plot_phase_portrait (self,t:float,n_X:int,n_Y:int, xlim:tuple, ylim:tuple,colors:list = None, fig = None, ax = None, plot_args = {},streamplot_args = {}):
        flag_ret_figs = False

        if fig == None and ax == None:
            fig,ax = plt.subplots()
            flag_ret_figs = True
        variable_x = self.variable_x
        variable_y = self.variable_y
        self.plot_streamplot(t,n_X,n_Y,xlim,ylim,fig,ax,**streamplot_args)

        if colors == None:
            self.plot_nullclines(variable_x,xlim,ylim,fig = fig, ax = ax,plot_args = plot_args)
            self.plot_nullclines(variable_y,xlim,ylim,fig = fig, ax = ax,plot_args = plot_args)
        else:
            self.plot_nullclines(variable_x,xlim,ylim,fig = fig, ax = ax,color = colors[0],plot_args = plot_args)
            self.plot_nullclines(variable_y,xlim,ylim,fig = fig, ax = ax,color = colors[1],plot_args = plot_args)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(self.variable_x)
        ax.set_ylabel(self.variable_y)
        if flag_ret_figs:
            return fig, ax
        
        
        
        
        
        
        
        