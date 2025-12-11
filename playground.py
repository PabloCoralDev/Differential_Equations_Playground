from dynamics_engine import DynamicsEngine

f = lambda x, t: -2*x**2 + 1

engine = DynamicsEngine(f, dt=0.1)
engine.get_response(initial_conditions=[0, 0], duration=10, return_as='plot')
