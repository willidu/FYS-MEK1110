import os
import numpy as np
import numpy.typing as npt
from typing import Optional, Tuple
import matplotlib.pyplot as plt         # type: ignore
from tqdm import trange                 # type: ignore

from potential import Lennard_Jones_Potential as LJP


class System:

    out_path = os.path.join(os.getcwd(), 'xyz_files')

    def __init__(
            self,
            r0: npt.ArrayLike,
            n: int,
            dim: int = 3,
            v0: Optional[npt.ArrayLike]=None, 
            L: Optional[float] = None,
            bound: bool = False,
            rc: Optional[float] = None,
            test: bool = False
        ) -> None:

        self.n = n
        self.dim = dim
        self.bound = bound
        self.rc = rc

        self.r0 = np.asarray(r0, dtype='float64')
        if v0 is None:
            self.v0 = np.zeros_like(self.r0)
        else:
            self.v0 = np.asarray(v0, dtype='float64')

        if type(L) is (int or float):
            self.L = L
        elif L is None and self.bound:
            self.L = 1.7*np.cbrt(self.n/4)

        def test_initials(self) -> None:
                if type(self.n) is not int:
                    raise IndexError(f'Number of particles must be integer, is {type(self.n)}')

                if len(self.r0) != self.n:
                    raise IndexError(
                        f'Incorrect length of positiolal vector for {self.n} atoms, length of r0 is {len(self.r0)}'
                    )

                if len(self.v0) != self.n:
                    raise IndexError(
                        f'Incorrect length of velocity vector for {self.n} atoms, length of v0 is {len(self.v0)}'
                    )

                if type(self.dim) is not int or self.dim<1 or self.dim>3:
                    raise IndexError(
                        f'Dimension must be integer between 1 and 3 (inclusive), is {type(self.dim)} with value {self.dim}'
                    )

                if self.dim >= 2:
                    for a in self.r0:
                        if len(a) != self.dim:
                            raise IndexError(
                                f'Incorrect dimensions for positional vectors, dimension is {len(a)}, should be {self.dim}'
                            )

                    for a in self.v0:
                        if len(a) != self.dim:
                            raise IndexError(
                                f'Incorrect dimensions for velocity vectors, dimension is {len(a)}, should be {self.dim}'
                            )

                if self.bound and (type(self.L) is (not int or  not float)):
                    raise TypeError(
                        f'Conflict in self.L and self.bound: {type(self.L) = }, {self.bound = }. self.L shound be int or float with boundaries turned on.'
                    )
                
                return

        if test:
            test_initials(self)

        print('System initiated successfully')

    def set_inital_velocities(self, v0: npt.ArrayLike=None, T: float=None) -> None:
        """ Sets initial velocity to v0 if v0 is explicit. 
        Else if temperature is zero, sets initial velocity to zero for all particles. 
        Else if T is nonzero, initial velocities get normally distributed with var=sqrt(T). 
        Else does nothing.
        """

        if v0 is not None:
            self.v0 = np.asarray(v0)

            if self.v0.shape != (self.n, self.dim):
                raise ValueError(
                    f'Incorrect dimension for initial velocities: should be ({self.n, self.dim}), is {v0.shape}'
                )

        elif T is not None and T != 0:
            self.v0 = np.random.normal(0, np.sqrt(T/119.7), (self.n, self.dim))

        elif T == 0:
            self.v0 = np.zeros((self.n, self.dim))

    def get_temperatures(self) -> np.ndarray:
        """ Returns average temperature for simulation if solve() has been done, esle raises AttributeError. """

        return 119.7/(3*self.n)*np.einsum('ijk,ijk->i', self.v, self.v)

    def calculate_distances(self, x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
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

    def calculate_acceleration(self, x: np.ndarray) -> Tuple[np.ndarray, float]:

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

    def solve(self, T: float, dt: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ Using the Velocity Verlet integration method """

        self.dt = dt
        numpoints = int(np.ceil(T/dt))
        t = np.linspace(0, T, numpoints)
        x = np.zeros((numpoints, self.n, self.dim))
        v = np.zeros_like(x)
        x[0] = self.r0; v[0] = self.v0
        
        a_, ep_ = self.calculate_acceleration(x[0])
        ep = np.zeros_like(t)
        ep[0] = ep_

        for i in trange(numpoints-1):
            x[i+1] = x[i] + v[i]*dt + 0.5*a_*dt**2
            a_2, ep_ = self.calculate_acceleration(x[i+1])
            v[i+1] = v[i] + 0.5*(a_ + a_2)*dt
            
            if self.bound:
                x[i+1] = x[i+1] - np.floor(x[i+1]/self.L)*self.L

            ep[i+1] = ep_
            a_ = a_2

        self.t = t; self.x = x; self.v = v
        self.ep = ep; self.ek = 0.5*np.einsum('ijk,ijk->i', self.v, self.v)

        print('Simulation completed')
        return t, x, v

    def write_xyz_file(self, filename: str) -> None:
        with open(os.path.join(System.out_path, filename), 'w') as file:
            for r in self.x:
                file.write(f'{len(r)} \n')
                file.write(f'type  x  y  z\n')
                for r_ in r:
                    file.write(f'Ar   {r_[0]}  {r_[1]}  {r_[2]}\n')
        print('Finished writing file',filename)

    def energy(self, show: bool=False) -> None:
        plt.plot(self.t, self.ep, label='Potential energy')

        plt.plot(self.t, self.ek, label='Kinetic energy')

        plt.plot(self.t, self.ek+self.ep, label='Total energy')
        
        plt.legend(ncol=3, loc='upper right')
        plt.xlabel('t*')
        plt.ylabel('Energy')
        plt.grid()
        
        if show:
            plt.show()

    def calculate_velocity_correlation(self):
        return np.sum(np.einsum('ijk,jk->ij', self.v, self.v0)/np.einsum('ij,ij->i',self.v0, self.v0), axis=1)/self.n

    def diffusion_coefficient(self):
        cutoff_index = int(3*self.dt)
        A = self.calculate_velocity_correlation()
        i = 1/3*np.trapz(A, self.t)
        return i
