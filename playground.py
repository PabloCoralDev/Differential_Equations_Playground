from dynamics_engine import DynamicsEngine


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


f = lambda t, y: -2*(1/(y**(2))) + 1 #& y -> y(t)

engine = DynamicsEngine(f, dt=0.1)

engine.get_response(initial_conditions=[0, 1], duration=10, return_as='plot')

#FOUND A SUPER SILLY ERROR: If my inital conditions are 0,0 then x=0 and 1/x^2 = 0...
#Must be able to handle singularities
