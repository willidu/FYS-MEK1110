import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
out_path = os.path.join(os.getcwd(), 'xyz_files')

class System:
    def __init__(self, r0, v0, n, dim, L=1, rc=None, bound=False, test=False):
        self.n = n
        self.dim = dim
        self.L = L
        self.bound = bound
        self.rc = rc

        def test_func():
            if len(r0) != n:
                raise IndexError(f'Incorrect length of positiolal vector for {n} atoms, length of r0 is {len(r0)}')

            if len(v0) != n:
                raise IndexError(f'Incorrect length of velocity vector for {n} atoms, length of v0 is {len(v0)}')

            if dim > 2:
                for a in r0:
                    if len(a) != dim:
                        raise IndexError(f'Incorrect dimensions for positional vectors, dimension is {len(a)}, should be {dim}')

                for a in v0:
                    if len(a) != dim:
                        raise IndexError(f'Incorrect dimensions for velocity vectors, dimension is {len(a)}, should be {dim}')
            return

        test_func() if test else None

        self.r0 = np.asarray(r0)
        self.v0 = np.asarray(v0)

        print('System initiated successfully')

    def a(self, x, t):
        a = np.zeros((self.n, self.dim))
        p = np.zeros(self.n)

        for i in range(self.n):
            r = x[np.arange(self.n)!=i] - x[i]

            # print(f'{i = }\n{r = }')

            if self.bound:
                """
                term = np.around(r/self.L)*self.L
                print(f'{np.linalg.norm(r-term, axis=1) = }')
                r -= np.around(r/self.L)*self.L
                print(f'{np.linalg.norm(r, axis=1) = }')
                exit()
                pass
                """
                r -= np.around(r/self.L)*self.L

            # print(f'{r = }')

            r_norm = np.linalg.norm(r, axis=1)
            # print(f'{r_norm = }') if (r_norm == 0).any() else None

            p[i] = 0.5 * np.sum(self.potential(r_norm))

            r_ = np.where(
                # r_norm < 3 and r_norm > 0, 
                np.logical_and(r_norm<3, r_norm>0.1),
                -1*(24*(2*(r_norm)**(-12)-(r_norm)**(-6))/(r_norm)**2), 
                0)

            a[i] = np.einsum('i, ij -> j', r_, r)
        return a, np.sum(p)

    def potential(self, r, sigma=1, epsilon=1):
        s6 = (sigma/r)**6
        s12 = s6 * s6
        if self.rc is None:
            return 4*epsilon*(s12-s6)
        else:
            return np.where(r < 3, 4*epsilon*(s12-s6) - 4*epsilon*((sigma/self.rc)**12 - (sigma/self.rc)**6), 0)

    def solve(self, T, dt):
        # Using the Velocity Verlet integration method
        numpoints = int(np.ceil(T/dt))
        t = np.linspace(0, T, numpoints)
        x = np.zeros((numpoints, self.n, self.dim))
        v = np.zeros_like(x)
        x[0] = self.r0; v[0] = self.v0
        
        ep = np.zeros_like(t)
        ek = np.zeros_like(t)

        a_, ep_ = self.a(x[0], t[0])
        ep[0] = ep_
        ek[0] = 0.5*np.sum(v[0]**2)

        for i in trange(numpoints-1):
            x[i+1] = x[i] + v[i]*dt + 0.5*a_*dt**2
            a_2, ep_ = self.a(x[i+1], t[i+1])
            v[i+1] = v[i] + 0.5*(a_ + a_2)*dt
            
            if self.bound:
                x[i+1] = x[i+1] - np.floor(x[i+1]/self.L)*self.L

            ep[i+1] = ep_
            ek[i+1] = 0.5*np.sum(v[i+1]**2)

            a_ = a_2

        print('Simulation completed')
        self.t = t; self.x = x; self.v = v; self.ep = ep; self.ek = ek
        return t, x, v

    def write_xyz_file(self, filename):
        with open(os.path.join(out_path, filename), 'w') as file:
            for r in self.x:
                file.write(f'{len(r)} \n')
                file.write(f'type  x  y  z\n')
                for r_ in r:
                    file.write(f'Ar   {r_[0]}  {r_[1]}  {r_[2]}\n')
        print('Finished writing file',filename)

    def energy(self, show=False):
        plt.plot(self.t, self.ep, label='Potential energy')

        plt.plot(self.t, self.ek, label='Kinetic energy')

        plt.plot(self.t, self.ek+self.ep, label='Total energy')
        
        plt.legend(ncol=3, loc='upper right')
        plt.xlabel('t*')
        plt.ylabel('Energy')
        plt.grid()
        plt.show() if show else None
