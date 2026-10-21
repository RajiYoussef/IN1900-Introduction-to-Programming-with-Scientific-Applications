from SEIR import *      # Imports from the file used to solve 1.2 to use the classes

# Oppgave a
class RegionInteraction(Region):
    def __init__(self, name, S_0, E2_0, lat, long):
        super().__init__(name, S_0, E2_0)
        self.lat = lat*np.pi/180
        self.long = long*np.pi/180

    def distance(self, other):
        lat_i = self.lat; long_i = self.long
        lat_j = other.lat; long_j = other.long
        R_earth = 64
        argument = np.sin(lat_i)*np.sin(lat_j) + np.cos(lat_i)*np.cos(lat_j)*np.cos(abs(long_i- long_j))
        if argument > 1 or argument < 0:    # Checks if argument to arccos is between 0 and 1
            return 0.0 #'The regions are the same
        else:
            distance = np.arccos(argument) * R_earth
            return distance                 # Returns in units of 10**4 meter


if __name__ == '__main__':
    innlandet = RegionInteraction('Innlandet',S_0=371385, E2_0=0, lat=60.7945,long=11.0680)
    oslo = RegionInteraction('Oslo',S_0=693494,E2_0=100, lat=59.9, long=10.8)
    print(oslo.distance(innlandet))

# Oppgave b
class ProblemInteraction(ProblemSEIR):
    def __init__(self, region, area_name, beta, r_ia = 0.1, r_e2=1.25,lmbda_1=0.33, lmbda_2=0.5, p_a=0.4, mu=0.2):
        # Uses the constructor from the ProblemSEIR class
        super().__init__(region, beta, r_ia, r_e2, lmbda_1, lmbda_2, p_a, mu)
        self.area_name = area_name

    def get_population(self):
        # Loops over the regions given in the region lists, and adds the poulations together
        s = 0
        for i in self.region:
            s += i.population
        return s

    def set_initial_condition(self):
        l = []
        # Loops through the regions and appends the intial conditions to initial_condition list
        for i in self.region:
            l.extend([i.S0, i.E1_0, i.E2_0, i.I0, i.Ia_0, i.R_0])
        self.initial_condition = l

    def __call__(self, u, t):
        lmbda_1, r_e2, r_ia, beta = self.lmbda_1, self.r_e2, self.r_ia, self.beta
        lmbda_2, p_a, mu = self.lmbda_2, self.p_a, self.mu
        n = len(self.region)
        # SEIR_list[i] = [S_i, E1_i, E2_i, I_i, Ia_i, R_i]:
        SEIR_list = [u[i:i+6] for i in range(0, len(u), 6)]
        # Create separate lists containing E2 and Ia values:
        E2_list = [u[i] for i in range(2, len(u), 6)]
        Ia_list = [u[i] for i in range(4, len(u), 6)]

        derivative = []
        for i in range(n): # County we are in
            S, E1, E2, I, Ia, R = SEIR_list[i]
            dS = 0
            for j in range(n): # The other county
                E2_other = E2_list[j]
                Ia_other = Ia_list[j]
            # Uses the forumals in the text to calulate the derivatives and returns them as a list
                dS += -Ia_other*np.exp(-self.region[i].distance(self.region[j]))/self.region[j].population*r_ia*beta(t)*S
                dS += -r_e2*beta(t)*S*E2_other*np.exp(-self.region[i].distance(self.region[j]))/self.region[j].population

            N = self.region[i].population
            dS += -beta(t)*S*I/N
            dI = lmbda_2*E2 - mu*I
            dE1 = -dS - lmbda_1*E1
            dE2 = lmbda_1*(1-p_a)*E1 - lmbda_2*E2
            dIa = lmbda_1*p_a*E1 - mu*Ia
            dR = mu*(I + Ia)

            derivative.extend([dS, dE1, dE2, dI, dIa, dR])
        return derivative

    def solution(self, u, t):
        n = len(t)
        n_reg = len(self.region)
        self.t = t
        self.S = np.zeros(n)
        self.E1 = np.zeros(n)
        self.E2 = np.zeros(n)
        self.I = np.zeros(n)
        self.Ia = np.zeros(n)
        self.R = np.zeros(n)
        SEIR_list = [u[:, i:i+6] for i in range(0, n_reg*6, 6)]
        for part, SEIR in zip(self.region, SEIR_list):
            # Loops through the regions, and list of SEIR values and adds together the different
            # SEIR values from each region
            part.set_SEIR_values(SEIR, t)
            self.S += part.S
            self.E1 += part.E1
            self.E2 += part.E2
            self.I += part.I
            self.Ia += part.Ia
            self.R += part.R

    def plot(self):
        # Plots the different SEIR values
        plt.title(f'{self.area_name}')
        plt.xlabel('Time(days)')
        plt.ylabel('Population')
        plt.plot(self.t, self.S, 'midnightblue', label = 'S')
        plt.plot(self.t, self.I, 'orange',label = 'I')
        plt.plot(self.t, self.Ia, 'g',label = 'Ia')
        plt.plot(self.t, self.R,'r', label = 'R')

if __name__ == '__main__':
    problem = ProblemInteraction([oslo, innlandet], 'Norway_east', beta=0.5)
    print(problem.get_population())
    problem.set_initial_condition()
    print(problem.initial_condition) #non-nested list of length 12
    u = problem.initial_condition
    print(problem(u,0)) #list of length 12. Check that values make sense
    solver = SolverSEIR(problem,T=100,dt=1.0)
    solver.solve()
    problem.plot()
    plt.legend()
    plt.show()

"""
Terminal > python SEIR_interaction.py
1.0100809386280782
1064979
[693494, 0, 100, 0, 0, 0, 371385, 0, 0, 0, 0, 0]
[-62.49098896472576, 62.49098896472576, -50.0, 50.0, 0.0, 0.0, -12.187832324277785, 12.187832324277785, 0.0, 0.0, 0.0, 0.0]
(shows a graph)
"""
