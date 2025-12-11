import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class Integrator:
    """
    Handles numerical integration of Circuit ODEs.
    Refactored from Combo_ODE_pp.
    """
    def __init__(self, circuit, parameters_values):
        self.circuit = circuit
        self.parameters_values = parameters_values
    
    def fun_to_integrate (self,Y,t):
        x,y = Y
        dxdt = self.circuit.xdot
        if len(self.parameters_values)>0:
            for p in self.circuit.parameters:
                dxdt = dxdt.subs(p,self.parameters_values[p])
        dxdt = dxdt.subs({self.circuit.x:x,self.circuit.y:y})
        dydt = self.circuit.ydot
        if len(self.parameters_values)>0:
            for p in self.circuit.parameters:
                dydt = dydt.subs(p,self.parameters_values[p])
        dydt = dydt.subs({self.circuit.x:x,self.circuit.y:y})
        return [dxdt,dydt]

    def integrate (self,Y0,t_final,dt,continue_to_decay = False,tol = 1e-4, time_cap = False, time_to_cap = 30):
        t_final = t_final
        dt = dt
        t = np.arange(0,t_final+dt,dt)
        ret = odeint(self.fun_to_integrate,Y0,t)
        xs = ret[:,0]
        ys = ret[:,1]
        if continue_to_decay:
            import time
            start_time = time.perf_counter()
            flag_cont = False
            while not flag_cont:
                if time_cap:
                    if time.perf_counter()-start_time<=time_to_cap:
                        if abs(xs[-1]-xs[-2])>tol or abs(ys[-1]-ys[-2])>tol:
                            t_final = 1.5*t_final
                            ret = self.integrate(Y0,t_final,dt)
                            t = ret['t'].values
                            xs = ret['x'].values
                            ys = ret['y'].values
                        else:
                            flag_cont = True
                    else:
                        break
                else:
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

    def integrate_with_noise(self, Y0, t_final, dt, sigmax, sigmay=0):
        t = np.arange(0, t_final + dt, dt)
        n_steps = len(t)
        traj = np.zeros((n_steps, 2))
        traj[0, :] = Y0

        for i in range(1, n_steps):
            # Compute deterministic update
            dxdt, dydt = self.fun_to_integrate(traj[i-1], t[i-1])

            # Additive white noise on x and y
            noise_x = np.random.normal(loc=0.0, scale=sigmax)
            noise_y = np.random.normal(loc=0.0, scale=sigmay)
            x_next = traj[i-1, 0] + dxdt * dt + noise_x * np.sqrt(dt)
            y_next = traj[i-1, 1] + dydt * dt + noise_y * np.sqrt(dt)

            traj[i, 0] = x_next
            traj[i, 1] = y_next

        df_ret = pd.DataFrame({'t': t, 'x': traj[:, 0], 'y': traj[:, 1]})
        return df_ret

    def integrate_and_plot (self,Y0,t_final,dt,fig = None, axs = None,continue_to_decay = False,
        time_cap = False, time_to_cap = 30):
            ret = self.integrate(Y0,t_final,dt)
            xs = ret[:,0]
            ys = ret[:,1]
            if continue_to_decay:
                import time
                start_time = time.perf_counter()
                flag_cont = False
                while not flag_cont:
                    if time_cap:
                        if time.perf_counter()-start_time<=time_to_cap:
                            if (np.argmax(xs)==len(xs)-1) or (np.argmax(ys)==len(ys)-1):
                                t_final = 1.5*t_final
                                ret = self.integrate(Y0,t_final,dt)
                                xs = ret[:,0]
                                ys = ret[:,1]
                            else:
                                flag_cont = True
                    else:
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
