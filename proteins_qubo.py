from dimod import ConstrainedQuadraticModel
cqm = ConstrainedQuadraticModel()
import numpy as np
from dimod import Integer
from dimod import Binary
E_HH = -2.3
E_HP = -1
E_PP = 0
import random

def get_adjacent_cells(dims, x, y):
    cells = []
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            if i != j:
                coord = [x+i, y+j]
                if coord[0] > -1 and coord[0] < dims[0] and coord[1] > -1 and coord[1] < dims[1]:
                    cells.append(coord)
    return cells


#find which amino acids are hydrophobic, record how many there are for the constraint
SEQ = 'ASN LEU TYR ILE GLN TRP LEU LYS ASP GLY GLY PRO SER SER GLY ARG PRO PRO PRO SER'
hydrophillic = ['ASP', 'ARG', 'LYS', 'GLU', 'ASN', 'SER', 'GLN']
SEQ = SEQ.split()
sigma = []
n = 0
for r in SEQ:
    if r in hydrophillic:
        sigma.append(0)
    else:
        sigma.append(1)
        n = n + 1
print(n, ' 1s')
size = (4, 5)
#create all variables in the QUBO problem
r = np.empty(size, dtype=object)
for x in range(size[0]):
    for y in range(size[1]):
        r[x, y] = Binary('r'+str(x)+str(y))
energy = []
ones = []
import random
for x in range(size[0]):
    for y in range(size[1]):
        coords = get_adjacent_cells(size, x, y)
        used = False
        i = 4*x + y
        for neighbor in coords:
            this_thing = r[x, y]
            neigh = r[neighbor[0], neighbor[1]]
            j = 4 * neighbor[0] + neighbor[1]
            if np.random.uniform() > 1.0/8.0:
                energy.append(-1 * (this_thing - neigh)**2 + .1 * ((this_thing - 1) + (neigh - 1))**2 - 2.3 * (this_thing + neigh)**2)
        ones.append(r[x, y])
cqm.set_objective(sum(energy))
cqm.add_constraint(sum(ones) == n, label=f'molecules')
print('Running sampler')

#sample proteins
from dwave.system import LeapHybridCQMSampler
sampler = LeapHybridCQMSampler()
sampleset = sampler.sample_cqm(cqm, time_limit=5, label="SDK Examples - Bin Packing")  
feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
def binary_grid(a):
    b = np.zeros(size)
    for coord, item in a:
        b[int(coord[1]), int(coord[2])] = item
    return b

#write the sanity check lattices to a file
if len(feasible_sampleset):      
    best_energy = feasible_sampleset.first
    print(feasible_sampleset.__dict__)
    print("{} feasible solutions of {}.".format(len(feasible_sampleset), len(sampleset)))
    soln_no = 0
    for soln in feasible_sampleset:
        with open("solution_no"+str(soln_no)+".txt", "w") as text_file:
            text_file.write(str(binary_grid(soln.items())))
        soln_no = soln_no + 1

