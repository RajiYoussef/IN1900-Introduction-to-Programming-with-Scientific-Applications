from ODEsolver import *
import matplotlib.pyplot as plt

# Oppgave a
def SEIR(u,t):
    ''' The SEIR function itself. Returns the parameters for the function. '''
    beta = 0.5; r_ia = 0.1; r_e2=1.25
    lmbda_1=0.33; lmbda_2=0.5; p_a=0.4; mu=0.2

    S, E1, E2, I, Ia, R = u
    N = sum(u)  # Calulates total population
    # Calulates the derivative
    dS = -beta*S*I/N - r_ia*beta*S*Ia/N - r_e2*beta*S*E2/N
    dE1 = beta*S*I/N + r_ia*beta*S*Ia/N + r_e2*beta*S*E2/N - lmbda_1*E1
    dE2 = lmbda_1*(1-p_a)*E1 - lmbda_2*E2
    dI = lmbda_2*E2 - mu*I
    dIa = lmbda_1*p_a*E1 - mu*Ia
    dR = mu*(I + Ia)
    # Returns a list with the calulated derivative values
    return [dS, dE1, dE2, dI, dIa, dR]

# Oppgave b
def solve(T, dt, S_0, E2_0):
    initial = [S_0, 0, E2_0, 0, 0, 0]       # Defines inital conditions
    points = int(T*dt + 1)                  # Start at 0 so needs to add a extra point
    timepoints = np.linspace(0, T, points)  # Calulates the time points for the ODE
    solution = RungeKutta4(SEIR)            # Initialises the ODEsolver
    solution.set_initial_condition(initial) # Sets the inital conditions for the ODEsolver
    u, t = solution.solve(timepoints)       # Computes the solution for the given time points
    return u, t

# Oppgave c
def plot(u,t):
    # Indexes the values we want to plot from the solution
    S = u[:,0]
    I = u[:,3]
    Ia = u[:,4]
    R = u[:,5]
    # Plots the values
    plt.plot(t, S, 'midnightblue', label = 'S')
    plt.plot(t, I, 'orange', label = 'I')
    plt.plot(t, Ia, 'g', label = 'Ia')
    plt.plot(t, R, 'r', label = 'R')
    plt.legend()
    plt.show()

# Oppgave d
# Sets inital conditions and plots
S_0 = 5e6; E2_0 = 100; T = 100; dt = 1.0
u, t = solve(T, dt, S_0, E2_0)
plot(u,t)
