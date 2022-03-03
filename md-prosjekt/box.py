import numpy as np

def unitcell(i, j, k, d):
    unit_cell = d * np.asarray(
            [
                [i, j, k], 
                [i, 0.5+j, 0.5+k],
                [0.5+i, j, 0.5+k],
                [0.5+i, 0.5+j, k]
            ]
        )
    return unit_cell

def lattice(n, L):
    d = L/n
    pos = np.zeros((4*n**3, 3))
    index = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                pos[index:index+4] = unitcell(i, j, k, d)
                index += 4
    return pos
