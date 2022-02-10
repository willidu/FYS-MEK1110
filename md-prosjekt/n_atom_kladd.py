import numpy as np
import matplotlib.pyplot as plt

class Atom:
    def __init__(self, r, n, dim):
        self.n = n
        self.dim = dim

        if len(r) != n:
            raise IndexError(f'Incorrect length of positiolal vector for {n} atoms, length of r is {len(r)}')
        for a in r:
            if len(a) != dim:
                raise IndexError(f'Incorrect dimensions for positional vectors, dimension is {len(a)}, should be {dim}')

        self.r = np.asarray(r)

    def a(self, t):
        # ri = self.r[0]
        for ri in self.r:
            for rj in self.r:
                if ri.all() == rj.all():
                    continue
                r_vec = np.array(ri-rj)
                r_norm = np.linalg.norm(r_vec)
                print(ri, rj, r_vec, r_norm)

        exit()

        term = 24*(2*r_norm**-12 - r_norm**-6)
        return term / (r_norm**2) * r_vec

    def Euler(self, N, T, t0, x0, v0):
        t = np.linspace(t0, T, N)
        x = np.zeros_like(t); x[0] = x0
        v = np.zeros_like(t); v[0] = v0
        a = np.zeros_like(t); a[0] = self.a
        dt = T/N

        for i in range(N+1):
            v[i+1] = v[i] + self.a * dt
            x[i+1] = x[i] + v[i] * dt

        return t, x, v, a

if __name__ == '__main__':
    pos_vec = [(3, -6), (-5, 2), (0, 5)]
    at1 = Atom(pos_vec, 3, 2)
    # at1 = Atom((3, -6))
    # at2 = Atom((-5, 2))
    r_ = at1.a(0)
    print(r_, len(r_), type(r_))

"""
class Atom:
    def __init__(self, r):
        self.r = np.asarray(r)

    def a(self, t):
        ri = self.r[0]
        rj = self.r[1]
        r_vec = np.array(ri-rj)
        r_norm = np.linalg.norm(r_vec)

        term = 24*(2*r_norm**-12 - r_norm**-6)
        vec = term / (r_norm**2) * r_vec
        avec = np.asarray([vec])
        return avec

    def Euler(self, T, dt):
        n = int(np.ceil(T/dt))
        t = np.linspace(0, T, n)
        x = np.zeros((n, len(self.r[0])))
        v = np.zeros_like(x)
        x[0] = self.r[0]

        for i in range(n-1):
            a_ = self.a(t[i])
            # print(type(a_), a_)
            v[i+1] = v[i] + self.a(t[i]) * dt
            x[i+1] = x[i] + v[i] * dt

        return t, x, v

pos_vec = ((0, 0), (0, 1.5))
atoms = Atom(pos_vec)
t, x, v = atoms.Euler(5, 0.01)
print(x)
"""