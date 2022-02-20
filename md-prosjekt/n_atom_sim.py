import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

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
        a = np.zeros((len(x), self.dim))
        p = np.zeros(len(x))

        for i in range(len(x)):
            r = x[np.arange(self.n)!=i] - x[i]
            
            if self.bound:
                r = r - round(r/self.L)*self.L

            r_norm = np.linalg.norm(r, axis=1)

            p[i] = np.sum(self.potential(r_norm))

            a[i] = np.sum(np.where(
                        r_norm < 3, 
                        -1*(24*(2*(r_norm)**(-12)-(r_norm)**(-6))*r/(r_norm)**2), 
                        0),
                    axis=0)
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
        ep = np.zeros_like(t)
        x[0] = self.r0; v[0] = self.v0

        a_, ep_ = self.a(x[0], t[0])
        ep[0] = ep_
        for i in range(numpoints-1):
            x[i+1] = x[i] + v[i]*dt + 0.5*a_*dt**2
            a_2, ep_ = self.a(x[i+1], t[i+1])
            ep[i+1] = ep_
            v[i+1] = v[i] + 0.5*(a_ + a_2)*dt
            a_ = a_2

        print('Simulation completed')
        self.t = t; self.x = x; self.v = v; self.ep = ep
        return t, x, v

    def write__xyz_file(self, filename, x):
        with open(filename, 'w') as file:
            for r in x:
                file.write(f'{len(r)} \n')
                file.write(f'type  x  y  z\n')
                for r_ in r:
                    file.write(f'Ar   {r_[0]}  {r_[1]}  {r_[2]}\n')
        print('Printing finished for file',filename)

    def energy(self, show=False):
        plt.plot(self.t, self.ep, label='Potential energy')

        # ek = np.sum(v**2) / 2
        ek = np.zeros_like(self.t)
        for i, v in enumerate(self.v):
            ek[i] = np.sum((np.dot(v_, v_) for v_ in v))

        plt.plot(self.t, ek, label='Kinetic energy')

        et = ek + self.ep
        plt.plot(self.t, et, label='Total energy')
        
        plt.legend(ncol=3, loc='upper right')
        plt.xlabel('t*')
        plt.ylabel('Energy')
        plt.grid()
        plt.show() if show else None
