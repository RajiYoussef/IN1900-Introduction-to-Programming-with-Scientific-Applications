import numpy as np
import matplotlib.pyplot as plt
from ODESolver import ExplicitMidpoint

# SEEIIR_model(u,t) is copied from the source and the name is changed to SEIR(u,t)
def SEIR(u,t):
    beta = 0.4
    r_ia = 0.1
    r_e2 = 1.25
    lmbda_1 = 0.33
    lmbda_2 = 0.5
    p_a = 0.4
    mu = 0.2

    S, E1, E2, I, Ia, R = u
    N = sum(u)
    dS  = -beta * S * I / N - r_ia * beta * S * Ia / N - r_e2 * beta * S * E2 / N
    dE1 = beta * S * I / N + r_ia * beta * S * Ia / N + r_e2 * beta * S * E2 / N - lmbda_1 *  E1
    dE2 = lmbda_1 * (1-p_a) * E1 - lmbda_2 * E2
    dI  = lmbda_2 * E2 - mu * I
    dIa = lmbda_1 * p_a*E1 - mu * Ia
    dR  = mu * (I + Ia)
    return [dS, dE1, dE2, dI, dIa, dR]

def test_SEIR():
    t = 0
    u = [1,1,1,1,1,1]
    expected = [-0.156666666666, -0.1733333333333, -0.302, 0.3, -0.068, 0.4]
    computed = SEIR(u,t)
    tol = 1E-10
    for i,j in zip(computed,expected):
        assert abs(i-j)<tol

    print(f"expected :{expected}")
    print(f"computed :{computed}")

test_SEIR()


def solve_SEIR(T,dt,S_0,E2_0):
    U0 = [S_0, 0, E2_0, 0, 0, 0]
    N = int(T / dt)
    time_points = np.linspace(0,T,N + 1)
    # i choose ExplicitMidpoint as a solver since the assignment did not mention any specific solver
    solver = ExplicitMidpoint(SEIR)
    solver.set_initial_condition(U0)
    u,t = solver.solve(time_points)
    return u,t


def plot_SEIR(u,t):
    plt.plot(t,u[:,0],label = "S(t)")
    plt.plot(t,u[:,3],label = "I(t)")
    plt.plot(t,u[:,4],label = "Ia(t)")
    plt.plot(t,u[:,5],label = "R(t)")
T = 150
dt = 1.0
S_0 = 5e6
E2_0 = 100
u,t = solve_SEIR(T,dt,S_0,E2_0)
plot_SEIR(u,t)
plt.title("Solved SEIR model with ExplicitMidpoint for categories,S,I,Ia,and,R")
plt.xlabel("Time (days)")
plt.ylabel("S,I,Ia,R")
plt.legend()
plt.grid()
plt.show()
"""
============================= test session starts =============================
collecting ... collected 1 item

seir_func.py::test_SEIR PASSED   [100%]expected :[-0.156666666666, -0.1733333333333, -0.302, 0.3, -0.068, 0.4]
computed : [-0.15666666666666668, -0.17333333333333334, -0.302, 0.3, -0.068, 0.4]


============================== 1 passed in 2.96s ==============================

Process finished with exit code 0

"""



