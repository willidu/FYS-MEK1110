import numpy as np
from n_atom_sim import System

def unitcell(i, j, k, d):
    unit_cell = d * np.asarray(
            [
                [i, j, k], 
                [i, 0.5+j, 0.5+k],
                [0.05+i, j, 0.5+k],
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

def task_3c():
    positions = lattice(n=3, L=20)

    with open('lattice_test.xyz', 'w') as file:
        file.write(f'{len(positions)} \n')
        file.write(f'type  x  y  z\n')
        for r_ in positions:
            file.write(f'Ar   {r_[0]}  {r_[1]}  {r_[2]}\n')

def task_3d():
    r0 = lattice(n=4, L=1.7*4)
    v0 = np.zeros_like(r0)
    s = System(r0, v0, 256, 3, True)
    print('System initiated')
    t, x, v = s.solve(5, 0.01)
    print('Simulation complete, writing to file.')
    s.write__xyz_file('task_3d.xyz', x)
    print('Energy calculations started')
    s.energy(show=True)
    print('Done!')

def task_3e():
    r0 = [1, 0, 0]
    v0 = [1, 0, 0]
    s = System(r0, v0, 1, 3, L=3, test=False)
    t, x, v = s.solve(15, 0.01)
    s.write__xyz_file('task_3e.xyz', x)

if __name__ == '__main__':
    # task_3c()
    # task_3d()
    task_3e()