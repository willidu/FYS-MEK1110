# kladd

import numpy as np
from n_atom_sim import *
from box import *
from potential import Lennard_Jones_Potential as LJP

r0 = lattice(4, 1.7*4)
v0 = np.zeros_like(r0)
v0[[i for i in range(1, 200, 7)]] = [1, .7, .1]
s = System(r0, v0, 256, 3, L=1.7*4, rc=3, bound=True, test=True)
s.solve(10, 0.005)
s.write_xyz_file('kladd.xyz')
s.energy(show=True)
print(s.ep)

# r0 = [[-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0]]
# r0 = np.asarray(r0, dtype='float64')
n = 1
L = 1.7 * n
r0 = lattice(n, L)
v0 = np.zeros_like(r0)
s = System(r0, v0, 4, 3, L=1.7, rc=3, bound=True, test=True)
s.solve(5, 0.01)
