# kladd

import numpy as np
from n_atom_sim import *
from box import *

s = MD(r0=lattice(n=2, L=1.7*2),
       n=32,
       bound=True,
       rc=3,
       test=True)
s.set_inital_velocities(T=300)
s.solve(5, 0.01)
s.write_xyz_file('kladd.xyz')
s.energy(True)
print(s.v[0])

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
