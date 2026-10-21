import matplotlib.pyplot as plt
class Region:
    '''Represents a geographical region.'''
    def __init__(self, name, S_0, E2_0):
        # Initalises the inital conditions for regions
        self.name = name
        self.S0 = S_0
        self.E2_0 = E2_0
        self.E1_0 = 0
        self.I0 = 0
        self.Ia_0 = 0
        self.R_0 = 0
        self.population = S_0 + E2_0

    def set_SEIR_values(self, u, t):
        # Assigns calulated values to the class
        self.S = u[:,0]; self.E1 = u[:,1]
        self.E2 = u[:,2]; self.I = u[:,3]
        self.Ia = u[:,4]; self.R = u[:,5]
        self.t = t

    def plot(self):
        # Plots the desired values
        S, I, Ia, R, t = self.S, self.I, self.Ia, self.R, self.t
        plt.xlabel('Time(days)')
        plt.ylabel('Population')
        plt.title(f'{self.name}')
        plt.plot(t, S, 'midnightblue', label = 'Suspectible')
        plt.plot(t, I, 'orange', label = 'Symptomatic infected')
        plt.plot(t, Ia, 'g', label = 'Asymptomatic infected')
        plt.plot(t, R, 'r',label = 'Recovered')