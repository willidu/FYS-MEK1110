# kladd

import numpy as np
from n_atom_sim import *
from box import *

d = np.zeros((50, int(5/0.01)))
for i in range(50):
       s = MD(r0=lattice(n=4, L=1.7*4),
              n=256,
              bound=True,
              rc=3,
              test=True)
       s.set_inital_velocities(T=180)
       t, x, v = s.solve(5, 0.01)
       d[i] = s.diffusion_coefficient()

plt.plot(t, np.average(d, axis=0), label='d')
plt.legend()
plt.show()