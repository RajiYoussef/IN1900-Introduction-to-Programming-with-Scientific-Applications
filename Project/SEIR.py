import numpy as np
import matplotlib.pyplot as plt
from ODESolver import *


class ProblemSEIR:
    def __init__(self, beta=0.4, r_ia=0.1, r_e2=1.25, \
                 lmbda_1=0.33, lmbda_2=0.5, p_a=0.4, mu=0.2):
        if isinstance(beta, (float, int)):
            self.beta = lambda t: beta
        elif callable(beta):
            self.beta = beta

        self.r_ia, self.r_e2, self.lmbda_1, self.lmbda_2, self.p_a, self.mu = r_ia, r_e2, lmbda_1, lmbda_2, p_a, mu

    def set_initial_condition(self, S_0, E2_0):
        self.S_0 = S_0
        self.E2_0 = E2_0
        initial_condition = [S_0, 0, E2_0, 0, 0, 0]
        self.initial_condition = initial_condition

    def get_population(self):

        return self.S_0 + self.E2_0

    def __call__(self, u, t):
        S, E1, E2, I, Ia, R = u
        N = sum(u)
        beta = self.beta(t)
        r_ia = self.r_ia
        r_e2 = self.r_e2
        lmbda_1 = self.lmbda_1
        lmbda_2 = self.lmbda_2
        p_a = self.p_a
        mu = self.mu

        dS = -beta * S * I / N - r_ia * beta * S * Ia / N - r_e2 * beta * S * E2 / N
        dE1 = beta * S * I / N + r_ia * beta * S * Ia / N + r_e2 * beta * S * E2 / N - lmbda_1 * E1
        dE2 = lmbda_1 * (1 - p_a) * E1 - lmbda_2 * E2
        dI = lmbda_2 * E2 - mu * I
        dIa = lmbda_1 * p_a * E1 - mu * Ia
        dR = mu * (I + Ia)
        return [dS, dE1, dE2, dI, dIa, dR]


def test_ProblemSEIR():
    S_0 = 1000  # random given values for the test
    E2_0 = 100
    # testing initial_condition
    test_problem = ProblemSEIR(beta=0.4, r_ia=0.1, r_e2=1.25, lmbda_1=0.33, lmbda_2=0.5, p_a=0.4, mu=0.2)
    tol = 1E-10
    test_problem.set_initial_condition(S_0, E2_0)
    computed = test_problem.initial_condition
    expected = [1000, 0, 100, 0, 0, 0]
    for i, j in zip(computed, expected):
        assert abs(i - j) < tol

    # testing get_population.Here i used a different method. the previous one gave me the following error :
    # for i, j in zip(computed, expected):
    #     E   TypeError: 'int' object is not iterable
    test_problem.get_population()
    computed = test_problem.get_population()
    expected = 1100
    tol = 1E-10
    success = abs(expected - computed) < tol
    msg = f"computed ={computed} {expected}=(expected)"
    assert success, msg

    # testing the call method
    t = 0
    u = [1, 1, 1, 1, 1, 1]
    expected = [-0.15666666666666668, -0.17333333333333334, -0.302, 0.3, -0.068, 0.4]
    computed = test_problem(u,t)
    tol = 1E-10
    for i,j in zip(computed,expected):
        assert abs(i - j) < tol


test_ProblemSEIR()


class SolverSEIR:
    def __init__(self,problem,T,dt):
        self.problem = problem
        self.T = T
        self.dt = dt

    def solve(self, method = ExplicitMidpoint):
        solver = method(self.problem)
        solver.set_initial_condition(self.problem.initial_condition)
        N = int(self.T / self.dt)
        t = np.linspace(0, self.T, N + 1)
        u, t = solver.solve(t)
        self.u, self.t = u,t

    def plot(self,states):
        states_dict = {'S': self.u[:, 0], 'I': self.u[:, 3], 'Ia': self.u[:, 4], 'R': self.u[:, 5]}
        for i in states_dict:
            plt.plot(self.t, states_dict[i], label=f"{i}(t)")



if __name__ == "__main__":
    test_ProblemSEIR()

S_0 = 5e6
E2_0 = 100
problem = ProblemSEIR(beta=0.4)
problem.set_initial_condition(S_0, E2_0)
solver = SolverSEIR(problem, T=150, dt=1.0)
solver.solve()
solver.plot(["S", "I", "Ia", "R"])
plt.title("Solved SEIR model with ExplicitMidpoint for categories,S,I,Ia,and,R")
plt.xlabel("Time(days)")
plt.ylabel("S,I,Ia,R")
plt.legend()
plt.grid()
plt.show()
"""
============================= test session starts =============================
collecting ... collected 2 items

SEIR.py::test_exact_numerical_solution <- ODESolver.py PASSED            [ 50%]
SEIR.py::test_ProblemSEIR PASSED                                         [100%]

============================== 2 passed in 2.50s ==============================

Process finished with exit code 0

"""