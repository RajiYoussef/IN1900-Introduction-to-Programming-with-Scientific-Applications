import matplotlib.pyplot as plt
from ODEsolver import *
from Region import Region 
from ProblemSEIR import ProblemSEIR as ps
from SolverSEIR import SolverSEIR as ss

# Oppgave a - Region
if __name__ == '__main__':
    nor = Region('Norway',S_0=5e6,E2_0=100)

    print(nor.name, nor.population)
    S_0, E1_0, E2_0 = nor.S0, nor.E1_0, nor.E2_0
    I_0, Ia_0, R_0 = nor.I0, nor.Ia_0, nor.R_0
    print(f'S_0 = {S_0}, E1_0 = {E1_0}, E2_0 = {E2_0}')
    print(f'I_0 = {I_0}, Ia_0 = {Ia_0}, R_0 = {R_0}')
    u = np.zeros((2,6))
    u[0,:] = [S_0, E1_0, E2_0, I_0, Ia_0, R_0]
    nor.set_SEIR_values(u,0)
    print(nor.S, nor.E1, nor.E2, nor.I, nor.Ia, nor.R)

# Oppgave b - ProblemSEIR
if __name__ == '__main__':  # Makes it possible to import from file
    problem = ps(nor,beta=0.5)
    problem.set_initial_condition()
    print(problem.initial_condition)
    print(problem.get_population())
    print(problem([1,1,1,1,1,1],0))

# Oppgave c - SolverSEIR
if __name__ == '__main__':  # Makes it possible to import from file
    solver = ss(problem,T=100,dt=1.0)
    solver.solve()
    nor.plot()
    plt.legend()
    plt.show()

