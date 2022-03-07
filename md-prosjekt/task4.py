from n_atom_sim import *
from box import *


def task_4_ii():
    r0 = lattice(n=3, L=1.7*3)
    
    s = System(r0, n=108, L=1.7*3, bound=True, rc=3, test=True)
    s.set_inital_velocities(T=300)
    t, x, v = s.solve(5, 0.001)
    s.write_xyz_file('4_a_ii.xyz')
    temp = s.get_temperature()
    print(temp[0])
    s.energy(True)

    plt.plot(t, temp, label='temp')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    task_4_ii()
