# Imports the classes used in 1.3
from SEIR_interaction import *

# Oppgave a
def read(filename):
    fylker = []
    with open(filename) as infile:
        for line in infile:
            # Loops through the text document, and stores the region from each line
            # as RegionInteraction class
            master = line.split(';')
            name = f'{master[1].strip()}'
            S_0 = float(master[2].strip())
            E2_0 = float(master[3].strip())
            lat = float(master[4].strip())
            long = float(master[5].strip())
            name = RegionInteraction(name, S_0, E2_0, lat, long)
            fylker.append(name)
    return fylker

# Oppgave b
def covid19_Norway(beta, filename, num_days, dt):
    # read file and create list of RegionInteraction instances
    list = read(filename)
    # create problem, an instance of ProblemInteraction
    problem = ProblemInteraction(list, 'Norway', beta)
    # create the solver, an instance of SolverSEIR
    solver = SolverSEIR(problem, num_days, dt)
    # call the method solve
    solver.solve()
    plt.figure(figsize=(9, 12)) # set figsize
    index = 1
    # for each part in problem’s attribute region:
    for i in problem.region:
        plt.subplot(4,3,index)
        # Calls plot method from current indexed region
        i.plot()
        index += 1

    plt.subplot(4,3, index)
    plt.subplots_adjust(hspace = 0.8, wspace=0.6)

    # Call plot method from problem
    problem.plot()
    # Makes a nice legend
    plt.subplot(431)
    plt.legend(bbox_to_anchor=(1., 1.02, 2., .102), loc='lower center', ncol= 1, mode="expand", borderaxespad=0.)
    # Shows the graph
    plt.show()

covid19_Norway(0.5, 'fylker.txt', 100, 1.0) # Calls the function

'''
0.2*0.05 = 0.01
The peak is around 500 000
500000*0.01 = 5000
There wouldn't be enough ventilators
'''

# Oppgave c
from datetime import *
def converter(R):
    return R/(1.25/0.5 + 0.1/0.2 + 1/0.2)

def piecewise_beta(t):
    # Calulates how many days from 15.2.2020 the given dates are
    t_1 = (date(2020, 3, 14) - date(2020, 2, 15)).days
    t_2 = (date(2020, 4, 20) - date(2020, 2, 15)).days
    t_3 = (date(2020, 5, 10) - date(2020, 2, 15)).days
    t_4 = (date(2020, 6, 30) - date(2020, 2, 15)).days
    t_5 = (date(2020, 7, 31) - date(2020, 2, 15)).days
    t_6 = (date(2020, 8, 30) - date(2020, 2, 15)).days
    t_7 = (date(2020, 9, 1) - date(2020, 2, 15)).days
    # Checks which of the intervals t is in, and returns a different R value based on that
    if t > 0 and t <= t_1:
        R = 4.0
    elif t > t_1 and t <= t_2:
        R = 0.5
    elif t > t_2 and t <= t_3:
        R = 0.4
    elif t > t_3 and t <= t_4:
        R = 0.8
    elif t > t_5 and t <= t_6:
        R = 0.9
    elif t > t_6 and t <= t_7:
        R = 1.0
    else:
        R = 1.1
    # Caluates the beta value from R, given the defintion in the text
    return R/(1.25/0.5 + 0.1/0.2 + 1/0.2)


time_passed = (date.today() - date(2020, 2, 15)).days
covid19_Norway(piecewise_beta, 'fylker.txt', time_passed, 1.0)
