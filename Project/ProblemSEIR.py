# Oppgave b
class ProblemSEIR:
    '''Defines the ODE model for a given region.'''
    def __init__(self, region, beta, r_ia = 0.1, r_e2=1.25, lmbda_1=0.33, lmbda_2=0.5, p_a=0.4, mu=0.2):
        # Checks if beta is constant or dependent on t
        if isinstance(beta, (float, int)):
            self.beta = lambda t: beta
        elif callable(beta):
            self.beta = beta

        self.region = region
        self.r_ia = r_ia
        self.r_e2 = r_e2
        self.lmbda_1 = lmbda_1
        self.lmbda_2 = lmbda_2
        self.p_a = p_a
        self.mu = mu

        # Sets the inital conditions by using a defined function
        self.set_initial_condition()

    def set_initial_condition(self):
        # Takes the intial values from the region in constructor and makes them
        # a part of ProblemSEIR as self.inital_condition
        reg = self.region
        self.initial_condition = [reg.S0, reg.E1_0, reg.E2_0, reg.I0, reg.Ia_0, reg.R_0]


    def get_population(self):
        # Gets the population from the region in the constructor
        return self.region.population

    def solution(self, u, t):
        return self.region.set_SEIR_values(u, t)

    def __call__(self, u, t):
        beta = self.beta; r_ia = self.r_ia; r_e2 = self.r_e2
        lmbda_1 = self.lmbda_1; lmbda_2 = self.lmbda_2; p_a = self.p_a; mu = self.mu

        S, E1, E2, I, Ia, R = u
        N = sum(u) # Calulates the population
        # Calulates the derivative values, and returns them as a list
        dS = -beta(t)*S*I/N - r_ia*beta(t)*S*Ia/N - r_e2*beta(t)*S*E2/N
        dE1 = beta(t)*S*I/N + r_ia*beta(t)*S*Ia/N + r_e2*beta(t)*S*E2/N - lmbda_1*E1
        dE2 = lmbda_1*(1-p_a)*E1 - lmbda_2*E2
        dI = lmbda_2*E2 - mu*I
        dIa = lmbda_1*p_a*E1 - mu*Ia
        dR = mu*(I + Ia)
        return [dS, dE1, dE2, dI, dIa, dR]