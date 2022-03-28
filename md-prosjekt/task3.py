from n_atom_sim import *
from box import *

def task_3a_iv():
    r = np.linspace(2, 4, 200)
    p = MD.LJP(r**2)
    p2 = MD.LJP(r**2, rc=3)

    plt.plot(r, p, label='Original potential')
    plt.plot(r, p2, label='shifted potential')
    plt.legend(loc='lower right')
    plt.xlim(2.5, 3.5)
    plt.xlabel('r*')
    plt.ylabel('Potential')
    plt.grid()
    plt.savefig(os.path.join(os.getcwd(), 'figures/3_a_iv.pdf'))
    plt.show()

def task_3b_i():
    s1 = MD(r0=[[0], [1.5]],
            n=2,
            dim=1,
            rc=3,
            test=True)
    t, x, v = s1.solve(5, 0.01)

    d = x[:,1] - x[:,0]
    plt.plot(t, d)
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.savefig(os.path.join(os.getcwd(), 'figures/3_b_i_15.pdf'))
    plt.show()

    s1.plot_energy()

    s2 = MD(r0=[[0], [0.95]],
            n=2,
            dim=1,
            rc=3,
            test=True)
    t, x, v = s2.solve(5, 0.01)

    d = x[:,1] - x[:,0]
    plt.plot(t, d)
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.savefig(os.path.join(os.getcwd(), 'figures/3_b_i_095.pdf'))
    plt.show()

    s2.plot_energy()

def task_3b_ii():
    s = MD(r0=[[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
           n=4,
           dim=3,
           rc=3,
           test=True)
    s.solve(5, 0.01)
    s.write_xyz_file('3_b_ii.xyz')

def task_3b_iv():
    s = MD(r0=[[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
           n=4,
           dim=3,
           rc=3,
           test=True)
    s.solve(5, 0.01)
    s.plot_energy()
    plt.savefig(os.path.join(os.getcwd(), 'figures/3_b_iv.pdf'))

def task_3b_v():
    s = MD(r0=[[1, 0.1, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
           n=4,
           dim=3,
           rc=3,
           test=True)
    s.solve(5, 0.01)
    s.write_xyz_file('3_b_v.xyz')
    s.plot_energy()
    plt.savefig(os.path.join(os.getcwd(), 'figures/3_b_v.pdf'))

def task_3c():
    positions = lattice(n=3, d=20/3)

    with open(os.path.join(MD.out_path, '3_c.xyz'), 'w') as file:
        file.write(f'{len(positions)} \n')
        file.write(f'type  x  y  z\n')
        for r_ in positions:
            file.write(f'Ar   {r_[0]}  {r_[1]}  {r_[2]}\n')

def task_3d():
    s = MD(r0=lattice(n=4),
           n=256,
           dim=3,
           rc=3,
           test=True)
    s.solve(5, 0.01)
    s.write_xyz_file('3_d.xyz')
    s.plot_energy()
    plt.savefig(os.path.join(os.getcwd(), 'figures/3_d.pdf'))

def task_3e():
    s = MD(r0=[1, 0, 0],
           n=1,
           dim=3,
           L=2,
           rc=3,
           bound=True)
    s.set_inital_velocities(v0=[1, 0, 0])
    s.solve(5, 0.01)
    s.write_xyz_file('3_e.xyz')

if __name__ == '__main__':
    task_3a_iv()
    task_3b_i()
    task_3b_ii()
    task_3b_iv()
    task_3b_v()
    task_3c()
    task_3d()
    task_3e()
