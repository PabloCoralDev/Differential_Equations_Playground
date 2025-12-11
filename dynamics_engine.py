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

    def get_response(self, initial_conditions=[0,1], duration=10, return_as='raw', tolerance=5) -> np.array:
        """
        READ ME!

        Simulate the response of the system for a specified time duration
        
        parameters: 
        - initial_coditions: float array, [t0, x0]. Defaults to [0,1]
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

        ##ERROR HANDLING ------------------------------------------

        """"
        try: #catch errors if function is not callable object 
            self.f(initial_conditions[0], initial_conditions[1])
        except Exception as e: #exception is a keyword
            raise ValueError("The function provided is not callable. Please provide a valid function f(x, t).") from e
        """

        ##mMAKE MY OWN RKF45 so I can finally understand how it works    
        ## Things to learn: f(t, y) returns a SLOPE at any point. Thus we can multiply by h or 1/2 h to get different solutions (slopes can be thought of as unit vectors, h is just length)

        ##INITIAL CONDITIONS ----------------------------------------
         
        #extract initial & final conditions: 
        t_0 = initial_conditions[0]
        t_end = duration
        y_0 = initial_conditions[1]

        t_pts = []
        y_pts = []
        s_factors = [1] #keep track to visualize later on. S here is INITIAL!
        

        num_pts = int((t_end - t_0) / self.dt) + 1 #+1 to include endpoint
        time_steps = np.linspace(t_0,t_end,num_pts) #need to include num_pts casted to int 
        #linspace is inclusive of start and endpoints. ALSO, 3d arg is NUMBER OF POINTS, not step size

        #calculate necessary coefficients
        y_k = y_0 #initial condition NECESSARY     

        for t_k in time_steps:

            ##LOCAL INITIAL CONDITIONS -----------------------------

            e = 1.0 #error starts at 100%
            e_cache = [e] #create a single array w errors
            s = s_factors[-1] #scaling factor starts at 1 ONLY on first step; Last accepted value of S is kept for the next run

            while True: 
                #for every step, run the calculation every time until the tolerance req is met

                #calculate the slope at the current y_k; I do NOT need t as an input!
                #time values will only dictate where to place the slope in the plane!
                #I can pass both t, y for nomenclature but t will have no effect.

                k1 = s * self.dt * self.f(t_k, y_k)
                #print(f'k1: {k1}')
                k2 = s * (self.dt) * self.f(t_k + (1/4)*self.dt, y_k + (1/4)*k1)
                k3 = s * (self.dt) * self.f(t_k + (3/8)*self.dt, y_k + (3/32)*k1 + (9/32)*k2)
                k4 = s * (self.dt) * self.f(t_k + (12/13)*self.dt, y_k + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3)
                k5 = s * (self.dt) * self.f(t_k + self.dt, y_k + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)
                k6 = s * (self.dt) * self.f(t_k + (1/2)*self.dt, y_k - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)

                #update y_new using weighted sum of ks. Create z_new for error estimation
                y_new = y_k + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5
                z_new = y_k + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6


                #recalculate error to ensure it's within tolerance
                #print(f'vals y:{y_new}, z: {z_new}')
                e_k = np.abs(z_new - y_new)
                e_cache.append(e_k)
                #print(f'\nerror {e_k}')

                #MUST CATCH ERRORS FROM ERROR!!!


                print(f'\ne_prev: {e_cache[0]} | \te_curr: {e_cache[1]} | \tdiff: {(np.abs(e_cache[0]-e_cache[1])/e_cache[1])} | \ttol: {tolerance}' )
                

                # OR if error converges is another case ; I should cache the error 
                # Apparently, relative error is STANDARD and goes alongisde the absolute error


                if e_k <= tolerance: 
                    y_k = y_new
                    t_pts.append(t_k)
                    y_pts.append(y_k)
                    break

                elif (np.abs(e_cache[0]-e_cache[1])/e_cache[1]) <= tolerance: #use the same tolerance for convergence, should be good enough
                    y_k = y_new
                    t_pts.append(t_k)
                    y_pts.append(y_k)
                    break

                #compare errors --> percent diff formula is (|curr-prev|/prev) -> denom is largest, which should be prev (prev is reference)
                #if no value is the reference, denominator becomes average (this is percent DIFFERENCE, not error) --> |curr-prev|/(|curr+prev|/2)
                else: 
                    s = 0.84 * (tolerance / e_k)**0.25 #formula for s given by Numerical Methods using MATLAB, 4th edition, Jaan Kiusalaas
                    s_factors.append(s)

                # if we did NOT exit the loop, I must adjust the cached errors before proceeding to the next run:
                e_cache.pop(0) #pop SHIFTS INDICES, which is what I want (leave only the current error.)
                # may have to do some sort of time-out in case tolerance is never met
                #if error is less than tolerance, loop will break naturally
                #if error is still greater than tolerance, step size needs to be recalculated
                  


            #program does NOT 'continue' if error is exceeded (step has to be adaptive from the beginning) 
            ## y_new is a 4th order estimate (up to k5), z_new is a 5th order estimate (up to k6). We can compare them to determine the next best step size

        if return_as == 'plot':

            plt.plot(t_pts, y_pts)
            plt.xlabel('Time')
            plt.ylabel('Response')
            plt.title('System Response Over Time')
            plt.grid()
            #plt.plot(t_pts, s_factors)
            plt.show()
            return None
        

        ## maybe one day I can implement this on vercel for a cool tool --> soon    
