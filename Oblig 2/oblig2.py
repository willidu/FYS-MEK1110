import numpy as np
import matplotlib.pyplot as plt

m = 0.1             # kg
L0 = 1              # m
k = 200             # N/m
g = 9.81            # m/s^2

theta0 = np.pi/6    # rad
r0 = L0*np.sin(theta0), -L0*np.cos(theta0)

T = 10
dt = 0.001
n = int(np.ceil(T/dt))

for k in (20, 2000):
    plt.clf()
    t = np.linspace(0, T, n)
    r = np.zeros((n, 2)); r[0] = r0
    v = np.zeros_like(r)
    a = np.zeros_like(r)

    for i in range(n-1):
        r_ = np.linalg.norm(r[i])
        a[i] = -k*(1-L0/r_)/m * r[i,0], -g-k*(1-L0/r_)/m * r[i, 1]
        v[i+1] = v[i] + a[i]*dt
        r[i+1] = r[i] + v[i+1]*dt

    plt.plot(r[:,0], r[:,1])
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.savefig(f'k{k}.pdf')
