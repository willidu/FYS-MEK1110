import numpy as np
import matplotlib.pyplot as plt

class System:
    def __init__(self, r0, v0, n, dim, test=False):
        self.n = n
        self.dim = dim

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

    def a(self, x, t):
        a = np.zeros((len(x), len(x), self.dim))
        for i in range(len(x)):
            for j in range(i+1, len(x)):
                r = x[i] - x[j]
                r_norm = np.linalg.norm(r)
                if r_norm <= 3:
                    a[i,j] = -(24*(2*(r_norm)**(-12) - (r_norm)**(-6)) * (r) / (r_norm)**2)
                    a[j,i] = -a[i,j]
                else:
                    a[i, j] = [0]*self.dim
                    a[j, i] = [0]*self.dim
        return np.sum(a, axis=0)

    def shifted_potential(self, r, cutoff=False, rc=None):
        # rc is cutoff point for distance r
        if cutoff and not rc is None:
            p = np.zeros_like(r)
            potential_at_cutoff = 4*(rc**(-12) - rc**(-6))
            for count, r_ in enumerate(r):
                if r_ < rc:
                    p[count] = 4*(r_**(-12) - r_**(-6)) - potential_at_cutoff
                else:
                    p[count] = 0
            return p
        else:
            return 4*(r**(-12) - r**(-6))

    def solve(self, T, dt):
        # Using the Velocity Verlet integration method
        numpoints = int(np.ceil(T/dt))
        t = np.linspace(0, T, numpoints)
        x = np.zeros((numpoints, self.n, self.dim))
        v = np.zeros_like(x)
        x[0] = self.r0; v[0] = self.v0

        a_ = self.a(x[0], t[0])
        for i in range(numpoints-1):
            x[i+1] = x[i] + v[i]*dt + 0.5*a_*dt**2
            a_2 = self.a(x[i+1], t[i+1])
            v[i+1] = v[i] + 0.5*(a_ + a_2)*dt
            a_ = a_2

        return t, x, v

    def write__xyz_file(self, filename, x):
        with open(filename, 'w') as file:
            for r in x:
                file.write(f'{len(r)} \n')
                file.write(f'type  x  y  z\n')
                for r_ in r:
                    file.write(f'Ar   {r_[0]}  {r_[1]}  {r_[2]}\n')
        print('printing done to',filename)

def task_3a_iv():
    r0 = [[0, 0], [1.5, 0]]
    v0 = np.zeros_like(r0)
    s1 = System(r0 , v0, 2, 2, test=True)
    t, x, v = s1.solve(5, 0.01)

    r = np.linspace(2, 4, 200)
    p = s1.shifted_potential(r)
    p2 = s1.shifted_potential(r, True, 3)

    plt.plot(r, p, label='Original potential')
    plt.plot(r, p2, label='shifted potential')
    plt.legend(loc='lower right')
    plt.xlim(2.5, 3.5)
    plt.xlabel('r*')
    plt.ylabel('Potential')
    plt.grid()
    plt.show()

def task_3b_i():
    r0 = [[0], [1.5]]
    v0 = np.zeros_like(r0)
    s1 = System(r0 , v0, 2, 1, test=True)
    t, x, v = s1.solve(5, 0.01)

    d = x[:,1] - x[:,0]
    plt.plot(t, d, label='$r_j - r_i$')
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.legend()
    plt.show()

    r0 = [[0], [0.95]]
    v0 = np.zeros_like(r0)
    s2 = System(r0 , v0, 2, 1, test=True)
    t, x, v = s2.solve(5, 0.01)

    d = x[:,1] - x[:,0]
    plt.plot(t, d, label='$r_j - r_i$')
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.legend()
    plt.show()

def task_3b_ii():
    r0 = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    v0 = np.zeros_like(r0)
    s1 = System(r0, v0, 4, 3, test=True)

    t, x, v = s1.solve(5, 0.001)
    s1.write__xyz_file('4atoms.xyz', x)

if __name__ == '__main__':
    pass