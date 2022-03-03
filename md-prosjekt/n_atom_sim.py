import os
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from potential import Lennard_Jones_Potential as LJP
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


class System:
    out_path = os.path.join(os.getcwd(), 'xyz_files')

    def __init__(self, r0, v0, n, dim, L=1, rc=None, bound=False, test=False):
        self.n = n
        self.dim = dim
        self.L = L
        self.bound = bound
        self.rc = rc

        def test_initials():
                if len(r0) != n:
                    raise IndexError(
                        f'Incorrect length of positiolal vector for {n} atoms, length of r0 is {len(r0)}'
                        )

                if len(v0) != n:
                    raise IndexError(
                        f'Incorrect length of velocity vector for {n} atoms, length of v0 is {len(v0)}'
                        )

                if dim >= 2:
                    for a in r0:
                        if len(a) != dim:
                            raise IndexError(
                                f'Incorrect dimensions for positional vectors, dimension is {len(a)}, should be {dim}'
                                )

                    for a in v0:
                        if len(a) != dim:
                            raise IndexError(
                                f'Incorrect dimensions for velocity vectors, dimension is {len(a)}, should be {dim}'
                                )
                return

        if test:
            test_initials()

        self.r0 = np.asarray(r0, dtype='float64')
        self.v0 = np.asarray(v0, dtype='float64')

        print('System initiated successfully')

    def calculate_distances(self, x: np.ndarray) -> np.ndarray:
        """
        dist : 
            Matrix with distance between all atoms in all dimensions. Note the sign.
        r_norm : 
            One-dimensional array with all relative distances.
        """

        dr = np.zeros((self.n, self.n, self.dim))
        
        for i in range(self.n):
            r = x[np.arange(self.n)!=i] - x[i]

            if self.bound:
                r -= np.around(r/self.L)*self.L

            dr[i, np.arange(self.n)!=i] = r

        r_norm_squared = np.einsum('ijk,ijk->ij', dr, dr)

        return dr, r_norm_squared

    def calculate_acceleration(self, x):
        dr, r_norm_squared = self.calculate_distances(x)

        potential_energy = 0.5 * np.sum(
            LJP.potential(
                r_norm_squared[r_norm_squared!=0],
                rc=self.rc
            )
        )

        s6 = r_norm_squared**(-3)
        s12 = s6 * s6

        force = np.where(
            np.logical_and(r_norm_squared<9, r_norm_squared!=0),
            -1*(24*(2*s12-s6)/r_norm_squared), 
            0
        )

        a = np.einsum('ij,ijk->ik', force, dr)

        return a, potential_energy

    def solve(self, T, dt):
        # Using the Velocity Verlet integration method
        numpoints = int(np.ceil(T/dt))
        t = np.linspace(0, T, numpoints)
        x = np.zeros((numpoints, self.n, self.dim))
        v = np.zeros_like(x)
        x[0] = self.r0; v[0] = self.v0
        
        ep = np.zeros_like(t)
        ek = np.zeros_like(t)

        a_, ep_ = self.calculate_acceleration(x[0])

        ep[0] = ep_
        ek[0] = 0.5*np.sum(v[0]**2)

        for i in trange(numpoints-1):
            x[i+1] = x[i] + v[i]*dt + 0.5*a_*dt**2
            a_2, ep_ = self.calculate_acceleration(x[i+1])
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
        with open(os.path.join(System.out_path, filename), 'w') as file:
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
        
        if show:
            plt.show()
