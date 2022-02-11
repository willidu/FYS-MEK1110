from n_atom_sim import *
from box import *

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

def task_3b_iv():
    r0 = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    v0 = np.zeros_like(r0)
    s1 = System(r0, v0, 4, 3, test=True)
    s1.solve(5, 0.001)
    s1.energy(show=True)

def task_3b_v():
    r0 = [[1, 0.1, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    v0 = np.zeros_like(r0)
    s1 = System(r0, v0, 4, 3, test=True)
    t, x, v = s1.solve(5, 0.001)
    s1.write__xyz_file('task3b.xyz', x)
    s1.energy(show=True)

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
    task_3a_iv()
    task_3b_i()
    task_3b_ii()
    task_3b_iv()
    task_3b_v()
    task_3c()
    task_3d()
    task_3e()