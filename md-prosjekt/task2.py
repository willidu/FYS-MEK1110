from two_atom_sim import *

def task_2b_i_ii():
    # 2b.i.
    r0 = [0, 1.5]
    v0 = [0, 0]
    t, x, v = EulerCromer(5, 0.01, r0, v0)

    # 2b.ii.
    d = x[:,1] - x[:,0]
    plt.plot(t, d, label='$r_j - r_i$')
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.legend()
    plt.show()

def task_2b_iv():
    r0 = [0, .95]
    v0 = [0, 0]
    t, x, v = VelocityVerlet(5, 0.01, r0, v0)
    d = x[:,1] - x[:,0]
    plt.plot(t, d, label='$r_j - r_i$')
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.legend()
    plt.show()

def task_2c_i():
    t, x, v = EulerCromer(5, 0.01, (0, 1.5), (0, 0))

    ep = np.asarray([U_marked(r) for r in x])
    plt.plot(t, ep, label='Potential energy')

    ek = np.asarray([np.dot(v_, v_) for v_ in v]) / 2
    plt.plot(t, ek, label='Kinetic energy')

    et = ek + ep
    plt.plot(t, et, label='Total energy')
    
    plt.legend(ncol=3, loc='upper right')
    plt.xlabel('t*')
    plt.ylabel('Energy')
    plt.ylim(-1.1, 1.)
    plt.grid()
    plt.show()

    t, x, v = EulerCromer(5, 0.01, (0, .95), (0, 0))

    ep = np.asarray([U_marked(r) for r in x])
    plt.plot(t, ep, label='Potential energy')

    ek = 0.5*v[:,0]**2 + 0.5*v[:,1]**2
    plt.plot(t, ek, label='Kinetic energy')

    et = ek + ep
    plt.plot(t, et, label='Total energy')
    
    plt.legend(ncol=3, loc='upper right')
    plt.xlabel('t*')
    plt.ylabel('Energy')
    plt.grid()
    plt.show()

def task_2c_iv():
    for dist in (1.5, .95):
        for solver in (Euler, EulerCromer, VelocityVerlet):
            t, x, v = solver(5, 0.01, (0, dist), (0, 0))
            ep = np.asarray([U_marked(r) for r in x])
            ek = 0.5*v[:,0]**2 + 0.5*v[:,1]**2
            et = ek + ep
            plt.plot(t, et, label=solver.__name__)
            plt.legend()
            plt.xlabel('t*')
            plt.ylabel('Energy')
            plt.grid()
        plt.show()

def task_2c_v():
    for dt in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
        print(dt)
        for solver in (EulerCromer, VelocityVerlet):
            t, x, v = solver(5, dt, (0, 1.5), (0, 0))
            plt.plot(t, x[:,0], label=f'xi, {solver.__name__}')
            plt.plot(t, x[:,1], label=f'xj, {solver.__name__}')
        
        plt.grid()
        plt.legend()
        plt.show()

        for solver in (EulerCromer, VelocityVerlet):
            t, x, v = solver(5, dt, (0, 1.5), (0, 0))
            ep = np.asarray([U_marked(r) for r in x])
            ek = 0.5*v[:,0]**2 + 0.5*v[:,1]**2
            et = ek + ep
            plt.plot(t, et, label=solver.__name__)
            plt.legend()
            plt.grid()
        plt.show()

def task_2d_i():
    t, x, v = VelocityVerlet(5, 0.01, (0, 1.5), (0, 0))
    with open('data1_5.xyz', 'w') as file:
        for r in x:
            file.write(f'{len(r)} \n')
            file.write(f'type  x  y  z\n')
            file.write(f'Ar   0  {r[0]}  0\n')
            file.write(f'Ar   0  {r[1]}  0\n')

    t, x, v = VelocityVerlet(5, 0.01, (0, .95), (0, 0))
    with open('data0_95.xyz', 'w') as file:
        for r in x:
            file.write(f'{len(r)} \n')
            file.write(f'type  x  y  z\n')
            file.write(f'Ar   0  {r[0]}  0\n')
            file.write(f'Ar   0  {r[1]}  0\n')


if __name__ == '__main__':
    task_2b_i_ii()
    task_2b_iv()
    task_2c_i()
    task_2c_iv()
    task_2c_v()
    task_2d_i()