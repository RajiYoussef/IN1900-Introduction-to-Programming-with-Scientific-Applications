from ODEsolver import *
class SolverSEIR:
    ''' Solves the SEIR system of ODEs for a given region.'''
    def __init__(self, problem, T, dt):
        # Initalises the values needed to solve the ODE
        self.problem = problem
        self.T = T
        self.dt = dt
        self.total_population = problem.get_population()

    def solve(self, method=RungeKutta4):
        # Solves the ODE and puts the values back in the ProblemSEIR class
        solver = method(self.problem)
        solver.set_initial_condition(self.problem.initial_condition)
        points = int(self.T*self.dt + 1)   # Includes the first value, so need to add one
        t = np.linspace(0, self.T, points)
        u, t = solver.solve(t)
        self.problem.solution(u, t)