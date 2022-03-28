from n_atom_sim import *
from box import *


def task_4_a_ii():
    s = MD(r0=lattice(n=3),
           n=108,
           bound=True,
           rc=3,
           test=True)
    s.set_inital_velocities(T=300)
    t, x, v = s.solve(5, 0.01)
    s.write_xyz_file('4_a_ii.xyz')
    temp = s.get_temperatures()

    plt.plot(t, temp)
    plt.xlabel('t*')
    plt.ylabel('Temperature [K]')
    plt.savefig(os.path.join(os.getcwd(), 'figures/4_a_ii.pdf'))
    plt.show()

def task_4_a_iii():
    dt = 0.01
    cutoff_index = int(0.5/dt)
    t_ = []
    for i in range(10):
        s = MD(r0=lattice(n=3),
               n=108,
               bound=True,
               rc=3,
               test=True)
        s.set_inital_velocities(T=180)
        t, x, v = s.solve(5, dt)
        temp = s.get_temperatures()
        temp = temp[cutoff_index:-1]
        t_.append(np.mean(temp))
    mean = np.mean(t_)
    print(mean)
    # 180 -> avg equlibium temp = 95.5

def task_4_b_ii():
    s = MD(r0=lattice(n=4),
           n=256,
           bound=True,
           rc=3,
           test=True)
    s.set_inital_velocities(T=180)
    t, x, v = s.solve(3, 0.01)

    A = s.vac()
    plt.plot(t, A)
    plt.xlabel('t*')
    plt.ylabel('Velocity autocorrelation')
    plt.savefig(os.path.join(os.getcwd(), 'figures/4_b_ii.pdf'))
    plt.show()

def task_4_b_iii():
    s1 = MD(r0=lattice(n=4),
            n=256,
            bound=True,
            rc=3,
            test=True)
    s1.set_inital_velocities(T=180)
    t, x, v = s1.solve(15, 0.01)

    s2 = MD(r0=x[-1],
            n=256,
            bound=True,
            rc=3,
            test=True)
    s2.set_inital_velocities(v0=v[-1])
    t, x, v = s2.solve(3, 0.01)
    
    A = s2.vac()
    plt.plot(t, A)
    plt.xlabel('t*')
    plt.ylabel('Velocity autocorrelation')
    plt.savefig(os.path.join(os.getcwd(), 'figures/4_b_iii.pdf'))
    plt.show()

def task_4_b_v():
    s1 = MD(r0=lattice(n=4),
            n=256,
            bound=True,
            rc=3,
            test=True)
    s1.set_inital_velocities(T=180)
    t, x, v = s1.solve(15, 0.01)

    s2 = MD(r0=x[-1],
            n=256,
            bound=True,
            rc=3,
            test=True)
    s2.set_inital_velocities(v0=v[-1])
    t, x, v = s2.solve(3, 0.01)
    
    D = s2.diffusion_coefficient()
    plt.plot(t, D)
    plt.xlabel('t*')
    plt.ylabel('Diffusion coefficient')
    plt.savefig(os.path.join(os.getcwd(), 'figures/4_b_v.pdf'))
    plt.show()
    
def task_4_c_ii():
    s = MD(r0=lattice(n=6),
           n=864,
           bound=True,
           rc=3,
           test=True)
    s.set_inital_velocities(T=180)
    s.solve(10, 0.01)

    msd, t = s.msd()
    plt.plot(t, msd)
    plt.xlabel('t*')
    plt.ylabel('Mean square displacement')
    plt.savefig(os.path.join(os.getcwd(), 'figures/4_c_ii.pdf'))
    plt.show()

def task_4_d_i():
    s = MD(r0=lattice(n=7),
           n=1372,
           bound=True,
           rc=3,
           test=True)
    s.set_inital_velocities(T=180)
    t, x, v = s.solve(10, 0.01)

    rdf_array, bin_centres = s.rdf(200)

    plt.axhline(y=1, color='black')
    plt.plot(bin_centres, rdf_array)
    plt.xlabel('r*')
    plt.ylabel('Radial distribution')
    plt.savefig(os.path.join(os.getcwd(), 'figures/4_d_i.pdf'))
    plt.show()
    
if __name__ == '__main__':
    task_4_a_ii()
    task_4_a_iii()
    task_4_b_ii()
    task_4_b_iii()
    task_4_b_v()
    task_4_c_ii()
    task_4_d_i()
