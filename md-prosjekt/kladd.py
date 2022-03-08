# kladd

import numpy as np
from n_atom_sim import *
from box import *

s = MD(r0=lattice(n=4, L=1.7*4),
       n=256,
       bound=True,
       rc=3,
       test=True)
s.set_inital_velocities(T=180)
t, x, v = s.solve(4, 0.01)
# s.write_xyz_file('kladd.xyz')
s.energy(True)
A = s.calculate_velocity_correlation()
plt.plot(t, A)
plt.show()
D = s.diffusion_coefficient()
print(f'{D:.2e}')
