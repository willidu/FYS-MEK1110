# kladd

import numpy as np
a = np.zeros((10, 2))
# print(a, a.shape)

a1 = np.array([2, 2])
a2 = np.array([0, 1])
# print(a1-a2)

r = np.asarray([(1, 1), (0, 0), (3, 3), (3, 4)])

for counter1, r_ in enumerate(r):
    for counter2, r__ in enumerate(r):
        if np.array_equal(r_, r__):
            continue

        if counter1 == counter2:
            print(counter1, counter2)

        print(r_, r__, r__-r_, counter1, counter2)

for c1, ri in enumerate(r):
    for c2, rj in enumerate(r):
        if c1 == c2:
            continue
        elif c2 == len(r)-1-c1:
            print(c1, c2)
        else: 
            print(c1, ri, c2, rj)

from two_atom_sim import *

dt = 1e-6

# t, x, v = Euler(5, dt, (0, 1.5), (0, 0))

# ep = np.asarray([U_marked(r) for r in x])
# plt.plot(t, ep, label='Potential energy')

# ek = 0.5*v[:,0]**2 + 0.5*v[:,1]**2
# plt.plot(t, ek, label='Kinetic energy')

# et = ek + ep
# plt.plot(t, et, label='Total energy')

# plt.legend(ncol=3, loc='upper right')
# plt.xlabel('t*')
# plt.ylabel('Energy')
# plt.grid()
# plt.show()

for solver in (Euler, EulerCromer, VelocityVerlet):
    t, x, v = solver(5, dt, (0, 1.5), (0, 0))
    ep = np.asarray([U_marked(r) for r in x])
    ek = 0.5*v[:,0]**2 + 0.5*v[:,1]**2
    et = ek + ep
    plt.plot(t, et, label=solver.__name__)
    plt.legend()
    plt.grid()
plt.show()