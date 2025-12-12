import numpy as np
import matplotlib.pyplot as plt


class DynamicsEngine:
    """A simple dynamics engine to simulate standard form ODEs using RK4(5), self-programmed
    for refreshing knowledge on numerical integration methods. Purposefully using NO GenAI or claude agent to keep brain sharp.
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
        - return_as: 'raw', 'plot', 'vector_field', 'interpolated' (default: 'raw')
            - raw returns a 2D np array with x-points and t-points
            - plot plots the response using matplotlib. No return.
            - vector field returns a vector field of eqn (slope field with unit length)
            - interpolated returns an akima1D object for response interpolation
                - Akima object can be plotted using matplotlib directly due to its built-in __call__(x) method (will call interpolated output for every input x)
        - tolerance: float, error tolerance for adaptive step sizing (default: 1e-6)

        return: depents on return_as arg. Can be None.
        """

        # Use differential equation! Can be ANY standard form ODE: dx/dt = f(x, t)

        ##ERROR HANDLING ------------------------------------------

        #1. self.f must be callable
        #2. differential equation must not be singular for the given IC's    

        ##LOCAL FUNCTIONS -----------------------------------------

        def plot_response(s_factors, t_pts, y_pts):

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

            # Left plot: System Response
            ax1.plot(t_pts, y_pts, color='#2E86AB', linewidth=2, label='y(t)')
            ax1.set_xlabel('Time (s)', fontsize=11, fontweight='bold')
            ax1.set_ylabel('Response', fontsize=11, fontweight='bold')
            ax1.set_title('System Response Over Time', fontsize=12, fontweight='bold', pad=15)
            ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.8)
            ax1.legend(loc='best', framealpha=0.9)

            # Right plot: Adaptive Step Size Factors
            step_indices = np.arange(len(s_factors))
            ax2.plot(step_indices, s_factors, color='#A23B72', linewidth=2, marker='o',
                     markersize=3, alpha=0.8, label='Step Size Factor')
            ax2.axhline(y=1.0, color='#F18F01', linestyle='--', linewidth=1.5,
                        alpha=0.7, label='Baseline (s=1)')
            ax2.set_xlabel('Iteration', fontsize=11, fontweight='bold')
            ax2.set_ylabel('Scaling Factor (s)', fontsize=11, fontweight='bold')
            ax2.set_title('Adaptive Step Size Evolution', fontsize=12, fontweight='bold', pad=15)
            ax2.grid(True, alpha=0.3, linestyle='--', linewidth=0.8)
            ax2.legend(loc='best', framealpha=0.9)

            plt.suptitle('RKF45 Dynamics Engine Simulation', fontsize=14, fontweight='bold', y=1.00)
            plt.show()
            return None
        
        def plot_vector_field():
            """
            Plots vector field of the thing. Complicated now that I think about it...
            """
            #create a grid with x by x points
            #evaluate RKF4(5) at higher and higher IC's to find the value at every point on the field
            #visualize as arrows (?)
            pass

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

                #MUST CATCH ERRORS FROM INCORRECT INPUT

                #print(f'\ne_prev: {e_cache[0]} | \te_curr: {e_cache[1]} | \tdiff: {(np.abs(e_cache[0]-e_cache[1])/e_cache[1])} | \ttol: {tolerance}' )
            

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
                
                s = 0.84 * (tolerance / e_k)**0.25 #formula for s given by Numerical Methods using MATLAB, 4th edition, Jaan Kiusalaas
                s_factors.append(s)

                # if we did NOT exit the loop, I must adjust the cached errors before proceeding to the next run:
                e_cache.pop(0) #pop SHIFTS INDICES, which is what I want (leave only the current error.)

                #MUST MAKE A TIME-OUT FEATURE if tolerance is not met in max iters

                print(f'length: {len(s_factors)}')
                

            #program does NOT 'continue' if error is exceeded (step has to be adaptive from the beginning) 
            ## y_new is a 4th order estimate (up to k5), z_new is a 5th order estimate (up to k6). We can compare them to determine the next best step size

        if return_as == 'plot': plot_response(s_factors=s_factors, t_pts=t_pts, y_pts=y_pts)
        elif return_as == 'vector_field': plot_vector_field()
        

        ## maybe one day I can implement this on vercel for a cool tool --> soon    
