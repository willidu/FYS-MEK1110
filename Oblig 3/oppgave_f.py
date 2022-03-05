import numpy as np
import matplotlib.pyplot as plt
import sys

m = 5       # kg
k = 500     # N/m
l0 = 0.5    # m
h = 0.3     # m
g = 9.81    # m/s^2

dt = 0.01
N = int(np.ceil(10/dt))

t = np.linspace(0, 10, N)
x = np.zeros_like(t)
v = np.zeros_like(t)

x[0] = float(sys.argv[1])

try:
    mu = float(sys.argv[2])
except IndexError:
    mu = None

for i in range(N-1):
    if mu is not None:
        N = m*g + k*h*(1-l0/np.sqrt(x[i]**2 + h**2))
        f = -mu*N*np.sign(v[i])
        a = (-k*x[i]*(1-l0/np.sqrt(x[i]**2 + h**2)) + f)/m
    else:
        a = -k/m*x[i]*(1-l0/np.sqrt(x[i]**2 + h**2))
    v[i+1] = v[i] + a*dt
    x[i+1] = x[i] + v[i+1]*dt

plt.subplot(2, 1, 1)
plt.plot(t, x, label='x(t)')
plt.ylabel('x [m]')
plt.subplot(2, 1, 2)
plt.plot(t, v, label='v(t)')
plt.ylabel('v [m/s]')
plt.xlabel('t [s]')
plt.show()

ek = 0.5*m*v**2
ep = m*g*h + 0.5*k*(np.sqrt(x**2+h**2)-l0)**2
plt.plot(x, ek, label='Kinetic energy')
plt.plot(x, ep, label='Potential energy')
plt.plot(x, ek+ep, label='total energy')
plt.ylabel('Energy [J]')
plt.xlabel('x [m]')
plt.legend()
plt.show()
