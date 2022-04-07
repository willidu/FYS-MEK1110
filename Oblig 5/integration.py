import numpy as np
import matplotlib.pyplot as plt

R = 0.15    # m
h = 1.0     # m
k = 10e3    # N/m
mu = 0.30
g = 9.81    # m/s^2
m = 1.0     # kg
v0x = 3.0   # m/s
dt = 0.001  # s

n = int(np.ceil(1/dt))
t = np.linspace(0, 1, n)
r = np.zeros((n, 2))
v = np.zeros((n, 2))
w = np.zeros(n)

r[0,1] = h
v[0] = v0x, 0

for i in range(n-1):
    N = k*(R-r[i,1])**(3/2) if r[i,1]<R else 0

    a = np.asarray([-mu*N/m, N/m - g])
    v[i+1] = v[i] + a * dt
    r[i+1] = r[i] + v[i+1] * dt

    alpha = -3/2 * mu*N/(m*R)
    w[i+1] = w[i] + alpha * dt

plt.plot(r[:,0], r[:,1])
plt.axis('equal')
plt.xlabel('Horizontal position [m]')
plt.ylabel('Vertical position [m]')
plt.grid()
plt.xlim(0, max(r[:,0]))
plt.savefig('position.pdf')
plt.show()

fig, ax = plt.subplots(3)

ax[0].plot(t, v[:,0], label='$v_x$')
ax[0].set_ylabel('$v_x$ [m s$^{-1}$]')

ax[1].plot(t, v[:,1], label='$v_y$')
ax[1].set_ylabel('$v_y$ [m s$^{-1}$]')

ax[2].plot(t, w, label='$\omega$')
ax[2].set_xlabel('Time [s]')
ax[2].set_ylabel('$\omega$ [rad s$^{-1}$]')

for ax_ in ax:
    ax_.legend()
    ax_.set_xlim(0, 1)
    ax_.grid()

plt.savefig('velocities.pdf')
plt.show()
