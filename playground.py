from dynamics_engine import DynamicsEngine
import numpy as np


"""
y is ANOTHER LAMBDA FUNCTION:

-> y(t) = a separate equation
- I can or cannot know y(t) analytically which is fine, so long as I have real INITIAL CONDITIONS for every y^n(t) and how it couples with dy/dt
- And also, if I knew y(t) analytically why would I need to do a numerical solver...

Basically: 
1. Know the relationship that an unknown function, y(t) makes with its derivative dy/dt
2. Know that the function y(t) exists at some initial condition t=0 for example
3. Get the value of y(t) at that initial condition, y(0) = 1 for example (y(0) = 0 is invalid for a SOME diff eq's)
4. Now I have the IC's for my differential equation : f(y,t) = f(1, 0) and I can use this to numerically find y(t
)
-> f(y(t), t) = the differential equation whih DEPENDS on function values of y and outputs a rate 

"""

#throw a BUNCH of equations to my engine, works most of the time!

f = lambda t, y: -.2*y #& y -> y(t) = this one WORKS (basic decay DE)
g = lambda t, y: y*(1-(y/10)) #basic (LOGISTIC ODE)
k = lambda t, y: y*(y**2-y/100) #random eqn
sinusoidal = lambda t, y: np.sin(y) + np.cos(y)
harder_eqn = lambda t, y: -10*((y*t)-4*np.sin(t+1)) - 5*np.cos(t+5) #super interesting case
gauss_pulse = lambda t, y: -y + 500*np.exp(-200*(t-2)**2)

f_engine = DynamicsEngine(f, dt=0.001)
g_engine = DynamicsEngine(g, dt=.01)
k_engine = DynamicsEngine(k, dt=.0001)
sinusoidal_engine = DynamicsEngine(sinusoidal, dt=0.01)
harder_eqn_engine = DynamicsEngine(harder_eqn, dt=.001)  
gauss_pulse_engine = DynamicsEngine(gauss_pulse, dt=.01)

#f_engine.get_response(initial_conditions=[0, 25], duration=5, return_as='plot', tolerance=.8e-14)
#g_engine.get_response(initial_conditions=[0, 1], duration=20, return_as='plot')
#k_engine.get_response(initial_conditions=[10, 12], duration=5, return_as='plot')
#sinusoidal_engine.get_response(initial_conditions=[0, np.pi/2], duration=5, return_as='plot')
harder_eqn_engine.get_response(initial_conditions=[4, 1], duration=50, return_as='plot', tolerance=.1)
gauss_pulse_engine.get_response(initial_conditions=[0, 5], duration=15, return_as='plot', tolerance=.0000001)

#FOUND A SUPER SILLY ERROR: If my inital conditions are 0,0 then x=0 and 1/x^2 = 0...
#Must be able to handle singularities
