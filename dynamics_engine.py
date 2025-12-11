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

    def get_response(self, initial_conditions=[0,0], duration=10, return_as='raw') -> np.array:
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

        return: depents on return_as arg. Can be None.
        """

        # Use differential equation! Can be ANY standard form ODE: dx/dt = f(x, t)

        try: #catch errors if function is not callable object 
            lambda func: self.f(initial_conditions[0], initial_conditions[1])
        except Exception as e: #exception is a keyword
            raise ValueError("The function provided is not callable. Please provide a valid function f(x, t).") from e
        
        if return_as == 'plot':
            t_points = np.arange(initial_conditions[0], initial_conditions[0] + duration, self.dt)
            x_points = np.zeros_like(t_points)
            x_points[0] = initial_conditions[1]

            for i in range(1, len(t_points)):
                t = t_points[i-1]
                x = x_points[i-1]

                k1 = self.f(x, t)
                k2 = self.f(x + 0.5 * self.dt * k1, t + 0.5 * self.dt)
                k3 = self.f(x + 0.5 * self.dt * k2, t + 0.5 * self.dt)
                k4 = self.f(x + self.dt * k3, t + self.dt)

                x_points[i] = x + (self.dt / 6) * (k1 + 2*k2 + 2*k3 + k4)

            plt.plot(t_points, x_points)
            plt.xlabel('Time')
            plt.ylabel('Response')
            plt.title('System Response Over Time')
            plt.grid()
            plt.show()
            return None
