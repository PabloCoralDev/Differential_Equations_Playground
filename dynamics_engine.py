import numpy as np
import matplotlib.pyplot as plt


class DynamicsEngine:
    """A simple dynamics engine to simulate standard form ODEs using RK4, self-programmed
    for refreshing knowledge on numerical integration methods.
    """
    def __init__(self, f, dt):
        """
        Initialize the dynamics engine.

        Parameters:
        f : Callable. The function defining the ODEs, should take state and time as inputs f(x, t)
        dt : float. The time step for integration.
        """
        self.f = f
        self.dt = dt

    def get_response(self, initial_conditions=[0,0], duration=10, return_as='raw', tolerance=1e-6) -> np.array:
        """
        READ ME!

        Simulate the response of the system for a specified time duration
        
        parameters: 
        - initial_coditions: float array, [t0, x0]. Defaults to [0,0]
        - duration: integer, simulation duration in seconds. Defaults to 10s
        - return_as: 'raw', 'plot' or 'interpolated' (default: 'raw')
            - raw returns a 2D np array with x-points and t-points
            - plot plots the response using matplotlib. No return.
            - interpolated returns an akima1D object for response interpolation
                - Akima object can be plotted using matplotlib directly due to its built-in __call__(x) method (will call interpolated output for every input x) 
            - tolerance: float, error tolerance for adaptive step sizing (default: 1e-6)

        return: depents on return_as arg. Can be None.
        """

        # Use differential equation! Can be ANY standard form ODE: dx/dt = f(x, t)

        # catch errors better; no initial condition => RKF45 will not work
        try: #catch errors if function is not callable object 
            lambda func: self.f(initial_conditions[0], initial_conditions[1])
        except Exception as e: #exception is a keyword
            raise ValueError("The function provided is not callable. Please provide a valid function f(x, t).") from e
        
        #extract initial & final conditions: 
        t_0 = initial_conditions[0]
        t_end = duration
        y_0 = initial_conditions[1]
        
        ## Code my own RK4 or maybe even a more advanced plotting method.--> MAKE MY OWN RKF45 so I can finally understand how it works    
        ## Things to learn: f(t, y) returns a SLOPE at any point. Thus we can multiply by h or 1/2 h to get different solutions (slopes can be thought of as unit vectors, h is just length)

        #first get number of necessary steps to satisfy duration with dt (boundary contidion)
        #generate linearly spaced time points from 0 to duration with step dt

        #ONE DAY I could try to implement a 'run until convergence or divergence' method, but for now just do fixed duration
        time_steps = np.linspace(t_0,t_end,self.dt) #linspace is inclusive of start and endpoints

        #calculate necessary coefficients
        y_k = y_0 #initial condition NECESSARY     

        for t_k in time_steps:
            
            k1 = self.dt * self.f(t_k, y_k)
            k2 = self.dt * self.f(t_k + (1/4)*self.dt, y_k + (1/4)*k1)
            k3 = self.dt * self.f(t_k + (3/8)*self.dt, y_k + (3/32)*k1 + (9/32)*k2)
            k4 = self.dt * self.f(t_k + (12/13)*self.dt, y_k + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3)
            k5 = self.dt * self.f(t_k + self.dt, y_k + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)
            k6 = self.dt * self.f(t_k + (1/2)*self.dt, y_k - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)

            #update y_new using weighted sum of ks. Create z_new for error estimation
            y_new = y_k + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5
            z_new = y_k + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6

            #check error & ensure it's within tolerance
            err = np.abs(z_new - y_new)

            if err <= tolerance: 

            #program does NOT 'continue' if error is exceeded (step has to be adaptive from the beginning) 
            ## y_new is a 4th order estimate (up to k5), z_new is a 5th order estimate (up to k6). We can compare them to determine the next best step size

        if return_as == 'plot':

            plt.plot(t_points, x_points)
            plt.xlabel('Time')
            plt.ylabel('Response')
            plt.title('System Response Over Time')
            plt.grid()
            plt.show()
            return None
        

        ## maybe one day I can implement this on vercel for a cool tool --> soon    
