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

r0 = ((1, 2), (3, 3))
v0 = ((0, 0), (0, 0))
s1 = System(r0 , v0, 3, 2, test=True)
