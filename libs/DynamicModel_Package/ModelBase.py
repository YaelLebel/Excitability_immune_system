import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

class DynamicModel:
    def __init__(self):
        self.variables = {'t':{'derivative_function':lambda variables,parameters:1,'parameters':[]}}  # Dictionary to store variables and their functions
        
    def add_variable(self, name:str, derivative_function:"function", parameters:dict):
        """
        Add a variable to the model with its time derivative function and parameters.
        
        Parameters:
            name (str): Name of the variable.
            derivative_function (function): Function describing the time derivative of the variable.
            parameters (dict): Parameters required by the derivative function.
        """
        self.variables[name] = {
            'derivative_function': derivative_function,
            'parameters': parameters
        }
    
    def compute_derivative(self, variable_values:dict):
        """
        Compute the time derivatives of all variables based on their functions and current values.
        
        Parameters:
            variable_values (dict): Dictionary containing current values of variables.
        
        Returns:
            dict: Dictionary containing computed derivatives for each variable.
        """
        derivatives = {}
        for variable, info in self.variables.items():
            derivative_function = info['derivative_function']
            parameters = info['parameters']
            # Call derivative function with current variable values and parameters
            derivatives[variable] = derivative_function(variable_values, parameters)
        return derivatives

    def euler_integrate (self,variables_initial:dict,t_final:float,dt:float):
        """
        Parameters
        ----------
        variables_initial : dictionary
            initial values for integration of all variables.
        t_final : float
            final time for integration.
        dt : float
            timestep for integration.

        Returns
        -------
        result_df : pd.DataFrame
            dataframe that contains values of all variables in all timepoints.

        """
        if 't' not in variables_initial.keys():
            variables_initial['t'] = 0
        result_df = pd.DataFrame(columns = list(self.variables.keys()))
        variables_initial_df = pd.DataFrame(variables_initial,index = [0])
        result_df = pd.concat([result_df,variables_initial_df],ignore_index=True)
        t = result_df['t'].values[-1] 
        while t<t_final:
            variable_values_t = result_df.iloc[-1].to_dict().copy()
            derivatives = self.compute_derivative(variable_values_t)
            result_values = {}
            for variable in derivatives.keys():
                result_values[variable] = variable_values_t[variable] + derivatives[variable]*dt
            t = result_values['t']
            result_values_df = pd.DataFrame(result_values,index=[0])
            result_df = pd.concat([result_df,result_values_df],ignore_index=True)
        
        return result_df
        
    def euler_integrate_keep_positive (self,variables_initial:dict,t_final:float,dt:float, tol:float = 1e-14):
        """
        Parameters
        ----------
        variables_initial : dictionary
            initial values for integration of all variables.
        t_final : float
            final time for integration.
        dt : float
            timestep for integration.

        Returns
        -------
        result_df : pd.DataFrame
            dataframe that contains values of all variables in all timepoints.

        """
        
        result_df = pd.DataFrame(columns = list(self.variables.keys()))
        variables_initial_df = pd.DataFrame(variables_initial,index = [0])
        result_df = pd.concat([result_df,variables_initial_df],ignore_index=True)
        t = result_df['t'].values[-1] 
        while t<t_final:
            variable_values_t = result_df.iloc[-1].to_dict().copy()
            derivatives = self.compute_derivative(variable_values_t)
            result_values = {}
            for variable in derivatives.keys():
                result_values[variable] = variable_values_t[variable] + derivatives[variable]*dt
                if variable_values_t[variable] + derivatives[variable]*dt<tol:
                    result_values[variable] = 0
            t = result_values['t']
            result_values_df = pd.DataFrame(result_values,index=[0])
            #result_values_df[result_values_df < tol ] = 0
            result_df = pd.concat([result_df,result_values_df],ignore_index=True)
        
        return result_df    

    def integrate_and_plot (self,variables_initial:dict,t_final:float,dt:float,keep_positive:bool,variables_to_plot:list,colors:list = [],scale:str = 'linear',
                            fig = None, axs = None):
        """
        Parameters
        ----------
        variables_initial : dictionary
            initial values for integration of all variables.
        t_final : float
            final time for integration.
        dt : float
            timestep for integration.

        keep_positive : Boolean
            wether or not to keep variables positive (non-negative) during integration.
        variables_to_plot : list
            list of variable names to plot.
        colors : list, optional
            list of colors for each variable. must be at least at the length of variables_to_plot. The default is [].
        colors : string, optional
            scale of y axis: linear or log The default is linear.
        results_df : pd.DataFrame, optional
            dataframe that stores result of integration. preferably returned from euler_integrate or euler_integrate_keep_positive.

        Returns
        -------
        time_integration:pd.DataFrame
            result of integration
        fig,axs.

        """
        if keep_positive:
            time_integration = self.euler_integrate_keep_positive(variables_initial,t_final,dt)
        else:
            time_integration = self.euler_integrate(variables_initial,t_final,dt)
        if fig == None and axs == None:
            fig,axs = plt.subplots(nrows = len(variables_to_plot))
            for i in range(len(variables_to_plot)):
                ax = axs[i]
                v = variables_to_plot[i]
                if len(colors)>0:
                    color = colors[i]
                    ax.plot(time_integration['t'],time_integration[v],label = v,color = color)
                else:
                    ax.plot(time_integration['t'],time_integration[v],label = v)
                ax.set_title(v)
                if scale=='log':
                    ax.set_yscale('log')
            fig.tight_layout()
            return time_integration,fig,axs
        else:
            for i in range(len(variables_to_plot)):
                ax = axs[i]
                v = variables_to_plot[i]
                if len(colors)>0:
                    color = colors[i]
                    ax.plot(time_integration['t'],time_integration[v],label = v,color = color)
                else:
                    ax.plot(time_integration['t'],time_integration[v],label = v)
                ax.set_title(v)
                if scale=='log':
                    ax.set_yscale('log')
            fig.tight_layout()
            return time_integration

    
    def plot_track_on_phase_space(self,df:pd.DataFrame,variables_to_plot:list,fig=None, ax = None, color = None, colormap = None,
                                  label:str = None):
        """
        Parameters
        ----------
        df : pd.DataFrame
            dataframe of time integration; preferably result of euler_integrate or euler_itegrate_keep_positive.
        variables_to_plot : list
            list of 2 variables to plot on phase space.
        fig : plt.fig, optional
            fig to plot on. The default is None.
        ax : plt.ax, optional
            ax to plot on. The default is None.
        color : string, optional
            color of scatter plot. The default is None.
        label : str, optional
            label of the scatter. The default is None 
        Returns
        -------
        fig,ax
        """
        if fig == None and ax == None:
            fig,ax = plt.subplots()
            
        if color != None:
            if label!=None:
                ax.scatter(df[variables_to_plot[0]],df[variables_to_plot[1]],color = color,zorder = 2,label = label)
            else:
                ax.scatter(df[variables_to_plot[0]],df[variables_to_plot[1]],color = color,zorder = 2)
        elif color == None and colormap != None:
            t = np.arange(len(df[variables_to_plot[0]]))
            if label!=None:
                im1 = ax.scatter(df[variables_to_plot[0]],df[variables_to_plot[1]],c=t,cmap = colormap,zorder = 2,label = label)
            else:
                im1 = ax.scatter(df[variables_to_plot[0]],df[variables_to_plot[1]],c=t,cmap = colormap,zorder = 2)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im1, cax=cax, orientation='vertical')
        else:
            if label!=None:
                ax.scatter(df[variables_to_plot[0]],df[variables_to_plot[1]],zorder = 2,label = label)
            else:
                ax.scatter(df[variables_to_plot[0]],df[variables_to_plot[1]],zorder = 2)

        ax.set_xlabel(variables_to_plot[0])
        ax.set_ylabel(variables_to_plot[1])
        if label!=None:
            ax.legend()
        return fig,ax
            
        
if __name__ == "__main__":
    # Example usage:
    def derivative_function_x(variables, parameters):
        return parameters['a'] * variables['y']  # Example function for x_dot
    
    def derivative_function_y(variables, parameters):
        return -parameters['b'] * variables['x']  # Example function for y_dot
    
    # Create an instance of the DynamicModel class
    model = DynamicModel()
    
    # Add variables and their derivative functions to the model
    model.add_variable('x', derivative_function_x, {'a': 1})
    model.add_variable('y', derivative_function_y, {'b': 100})
    
    # Example usage: Compute derivatives based on current variable values
    variable_values = {'t':0,'x': 1, 'y': 2}
    computed_derivatives = model.compute_derivative(variable_values)
    print(computed_derivatives)  # Output: {'x': 2, 'y': -1}

    df = model.euler_integrate({'t':0,'x':1,'y':10}, 1, 0.1)
    print(df)