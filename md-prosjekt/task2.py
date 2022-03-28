from two_atom_sim import *
import os


def task_2b_i_ii():
    # 2b.i.
    r0 = [0, 1.5]
    v0 = [0, 0]
    t, x, v = EulerCromer(5, 0.01, r0, v0)

    # 2b.ii.
    d = x[:,1] - x[:,0]
    plt.plot(t, d)
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.savefig(os.path.join(os.getcwd(), 'figures/2_b_ii.pdf'))
    plt.show()

def task_2b_iv():
    r0 = [0, .95]
    v0 = [0, 0]
    t, x, v = EulerCromer(5, 0.01, r0, v0)
    d = x[:,1] - x[:,0]
    plt.plot(t, d)
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.savefig(os.path.join(os.getcwd(), 'figures/2_b_iv.pdf'))
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
    plt.savefig(os.path.join(os.getcwd(), 'figures/2_c_i_1.pdf'))
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
    plt.savefig(os.path.join(os.getcwd(), 'figures/2_c_i_2.pdf'))
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
        name = (str(dist)).replace('.','')
        plt.savefig(os.path.join(os.getcwd(), f'figures/2_c_iv_{name}.pdf'))
        plt.show()

def task_2c_v():
    for dt in (1e-3, 1e-4):
        print(dt)
        for solver in (Euler, EulerCromer, VelocityVerlet):
            t, x, v = solver(5, dt, (0, 1.5), (0, 0))
            d = x[:,1] - x[:,0]
            plt.plot(t, d, label=solver.__name__)
        
        plt.legend()
        plt.xlabel('t*')
        plt.ylabel('r*')
        plt.grid()
        plt.savefig(os.path.join(os.getcwd(), f'figures/2_c_v_{dt}_dist.pdf'))
        plt.show()

        for solver in (Euler, EulerCromer, VelocityVerlet):
            t, x, v = solver(5, dt, (0, 1.5), (0, 0))
            ep = np.asarray([U_marked(r) for r in x])
            ek = 0.5*v[:,0]**2 + 0.5*v[:,1]**2
            et = ek + ep
            plt.plot(t, et, label=solver.__name__)
        plt.legend()
        plt.xlabel('t*')
        plt.ylabel('Energy')
        plt.grid()
        plt.savefig(os.path.join(os.getcwd(), f'figures/2_c_v_{dt}_energy.pdf'))
        plt.show()

def task_2d_i():
    t, x, v = VelocityVerlet(5, 0.01, (0, 1.5), (0, 0))
    with open(os.path.join(os.getcwd(), 'xyz_files/2_d_i_1.xyz'), 'w') as file:
        for r in x:
            file.write(f'{len(r)} \n')
            file.write(f'type  x  y  z\n')
            file.write(f'Ar   0  {r[0]}  0\n')
            file.write(f'Ar   0  {r[1]}  0\n')

    t, x, v = VelocityVerlet(5, 0.01, (0, .95), (0, 0))
    with open(os.path.join(os.getcwd(), 'xyz_files/2_d_i_2.xyz'), 'w') as file:
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