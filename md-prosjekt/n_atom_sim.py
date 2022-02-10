import numpy as np
import matplotlib.pyplot as plt
import numba

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
                a[i,j] = -(24*(2*(r_norm)**(-12) - (r_norm)**(-6)) * (r) / (r_norm)**2)
                a[j,i] = -a[i,j]
        return np.sum(a, axis=0)

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

if __name__ == '__main__':
    r0 = [[0, 0], [1.5, 0]]
    v0 = np.zeros_like(r0)
    s1 = System(r0 , v0, 2, 2, test=True)
    t, x, v = s1.solve(5, 0.01)

    print(x[0])

    d = x[:,1] - x[:,0]
    d = d[:,0]
    plt.plot(t, d)
    plt.show()
