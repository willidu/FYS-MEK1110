# kladd

import numpy as np
from n_atom_sim import *
from box import *

# r0 = [[1, 0, 0], [1, 1, 0], [-1, 0, 0], [0, -1, 0]]
r0 = lattice(n=1, L=1.7)
v0 = np.zeros_like(r0)
s = System(r0, v0, 4, 3, L=1.7, rc=3, bound=True, test=True)
# s.write_xyz_file('kladd.xyz')
s.solve(5, 0.01)

"""
for i in range(1):
    for j in range(i+1,len(r0)):
        dr = r0[i] - r0[j]
        dr_ = dr
        if s.bound:
            dr = dr - np.round(dr/1.7)*1.7
            print(np.linalg.norm(dr))
"""