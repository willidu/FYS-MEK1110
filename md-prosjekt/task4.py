from n_atom_sim import *
from box import *


def task_4_a_ii():
    s = System(r0=lattice(n=3, L=1.7*3),
               n=108,
               bound=True,
               rc=3,
               test=True)
    s.set_inital_velocities(T=300)
    t, x, v = s.solve(5, 0.001)
    s.write_xyz_file('4_a_ii.xyz')
    temp = s.get_temperatures()

    plt.plot(t, temp, label='temp')
    plt.legend()
    plt.show()

def task_4_a_ii():
    dt = 0.001
    cutoff_index = int(0.5/dt)
    t_ = []
    for i in range(10):
        s = System(r0=lattice(n=3, L=1.7*3),
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
    # 174 -> avg equlibium temp = 89.44

def task_4_b_ii():
    s = System(r0=lattice(n=4, L=1.7*4),
               n=256,
               bound=True,
               rc=3,
               test=True)
    s.set_inital_velocities(T=180)
    t, x, v = s.solve(3, 0.001)
    A = s.calculate_velocity_correlation()
    plt.plot(t, A)
    plt.show()

def task_4_b_iii():
    s1 = System(r0=lattice(n=4, L=1.7*4),
                n=256,
                bound=True,
                rc=3,
                test=True)
    s1.set_inital_velocities(T=180)
    t, x, v = s1.solve(5, 0.001)

    equil_pos = x[-1]
    equil_vel = v[-1]

    s2 = System(r0=equil_pos,
                n=256,
                bound=True,
                rc=3,
                test=True)
    s2.set_inital_velocities(v0=equil_vel)
    t, x, v = s2.solve(3, 0.001)
    A = s2.calculate_velocity_correlation()
    plt.plot(t, A)
    plt.show()

def task_4_b_v():
    s1 = System(r0=lattice(n=4, L=1.7*4),
                n=256,
                bound=True,
                rc=3,
                test=True)
    s1.set_inital_velocities(T=180)
    t, x, v = s1.solve(5, 0.001)

    equil_pos = x[-1]
    equil_vel = v[-1]

    s2 = System(r0=equil_pos,
                n=256,
                bound=True,
                rc=3,
                test=True)
    s2.set_inital_velocities(v0=equil_vel)
    s2.solve(3, 0.001)
    
    D = s2.diffusion_coefficient()
    print(f'{D:.2e}')
    # Kjøring gir 3.48e-02
    

if __name__ == '__main__':
    # task_4_a_ii()
    # task_4_a_ii()
    # task_4_b_ii()
    # task_4_b_iii()
    task_4_b_v()
