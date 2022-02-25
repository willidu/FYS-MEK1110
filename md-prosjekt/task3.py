from n_atom_sim import *
from box import *

def task_3a_iv():
    r0 = [[0, 0], [1.5, 0]]
    v0 = np.zeros_like(r0)
    s1 = System(r0 , v0, 2, 2, test=True)
    s2 = System(r0 , v0, 2, 2, rc=3, test=True)

    r = np.linspace(2, 4, 200)
    p = s1.potential(r)
    p2 = s2.potential(r)

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
    s1 = System(r0 , v0, 2, 1, rc=3, test=True)
    t, x, v = s1.solve(5, 0.01)

    d = x[:,1] - x[:,0]
    plt.plot(t, d, label='$r_j - r_i$')
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.legend()
    plt.show()

    s1.energy(show=True)

    r0 = [[0], [0.95]]
    v0 = np.zeros_like(r0)
    s2 = System(r0 , v0, 2, 1, rc=3, test=True)
    t, x, v = s2.solve(5, 0.01)

    d = x[:,1] - x[:,0]
    plt.plot(t, d, label='$r_j - r_i$')
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.legend()
    plt.show()

    s2.energy(show=True)

def task_3b_ii():
    r0 = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    v0 = np.zeros_like(r0)
    s1 = System(r0, v0, 4, 3, rc=3, test=True)

    s1.solve(5, 0.01)
    s1.write_xyz_file('3b_ii.xyz')

def task_3b_iv():
    r0 = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    v0 = np.zeros_like(r0)
    s1 = System(r0, v0, 4, 3, rc=3, test=True)
    s1.solve(5, 0.01)
    s1.energy(show=True)

def task_3b_v():
    r0 = [[1, 0.1, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    v0 = np.zeros_like(r0)
    s1 = System(r0, v0, 4, 3, rc=3, test=True)
    s1.solve(5, 0.01)
    s1.write_xyz_file('3b_v.xyz')
    s1.energy(show=True)

def task_3c():
    positions = lattice(n=3, L=20)

    with open(os.path.join(out_path, '3_c.xyz'), 'w') as file:
        file.write(f'{len(positions)} \n')
        file.write(f'type  x  y  z\n')
        for r_ in positions:
            file.write(f'Ar   {r_[0]}  {r_[1]}  {r_[2]}\n')

def task_3d():
    r0 = lattice(n=4, L=1.7*4)
    v0 = np.zeros_like(r0)
    s = System(r0, v0, 256, 3, True)
    s.solve(5, 0.01)
    s.write_xyz_file('3d.xyz')
    s.energy(show=True)

def task_3e():
    r0 = [1, 0, 0]
    v0 = [1, 0, 0]
    s = System(r0, v0, 1, 3, L=2, rc=3, bound=True, test=False)
    s.solve(5, 0.01)
    s.write_xyz_file('3e.xyz')

if __name__ == '__main__':
    task_3a_iv()
    task_3b_i()
    task_3b_ii()
    task_3b_iv()
    task_3b_v()
    task_3c()
    task_3d()
    task_3e()