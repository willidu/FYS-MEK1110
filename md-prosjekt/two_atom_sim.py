import numpy as np
import matplotlib.pyplot as plt

def a(r, t):
    ri = r[0]
    rj = r[1]
    r_vec = ri-rj
    r_norm = np.linalg.norm(r_vec)

    term = 24*(2*r_norm**-12 - r_norm**-6)
    vec = term / (r_norm**2) * r_vec
    return np.asarray([vec, -vec])

def Euler(T, dt, r0, v0):
    n = int(np.ceil(T/dt))
    t = np.linspace(0, T, n)
    x = np.zeros((n, 2))
    v = np.zeros_like(x)
    x[0] = r0; v[0] = v0

    for i in range(n-1):
        v[i+1] = v[i] + a(x[i], t[i]) * dt
        x[i+1] = x[i] + v[i] * dt

    return t, x, v

def EulerCromer(T, dt, r0, v0):
    n = int(np.ceil(T/dt))
    t = np.linspace(0, T, n)
    x = np.zeros((n, 2))
    v = np.zeros_like(x)
    x[0] = r0; v[0] = v0

    for i in range(n-1):
        v[i+1] = v[i] + a(x[i], t[i]) * dt
        x[i+1] = x[i] + v[i+1] * dt

    return t, x, v

def VelocityVerlet(T, dt, r0, v0):
    n = int(np.ceil(T/dt))
    t = np.linspace(0, T, n)
    x = np.zeros((n, 2))
    v = np.zeros_like(x)
    x[0] = r0; v[0] = v0

    a_ = a(x[0], t[0])
    for i in range(n-1):
        x[i+1] = x[i] + v[i]*dt + 0.5*a_*dt**2
        a_2 = a(x[i+1], t[i+1])
        v[i+1] = v[i] + 0.5*(a_ + a_2)*dt
        a_ = a_2

    return t, x, v

def U_marked(r):
    ri = r[0]
    rj = r[1]
    r_norm = abs(np.linalg.norm(rj-ri))
    return 4*(r_norm**(-12) - r_norm**(-6))

if __name__ == '__main__':
    # 2b.i.
    r0 = (0, 1.5)
    v0 = (0, 0)
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

    # 2.b.iv.
    r0 = (0, .95)
    v0 = (0, 0)
    t, x, v = VelocityVerlet(5, 0.01, r0, v0)
    d = x[:,1] - x[:,0]
    plt.plot(t, d, label='$r_j - r_i$')
    plt.xlabel('t*')
    plt.ylabel('r*')
    plt.ylim(0.8)
    plt.grid()
    plt.legend()
    plt.show()

    # 2.c.i.
    t, x, v = EulerCromer(5, 0.01, (0, 1.5), (0, 0))

    ep = np.asarray([U_marked(r) for r in x])
    plt.plot(t, ep, label='Potential energy')

    ek = 0.5*v[:,0]**2 + 0.5*v[:,1]**2
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

    # 2.c.iv.
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

    # 2.c.v.
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

    # 2.di.
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
