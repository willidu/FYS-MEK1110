import numpy as np
import matplotlib.pyplot as plt
import sys


"""
Kjøres i terminal.                   x   mu(valgfri)
Uten friksjon: > python oppgave_f.py 0.6
Med friksjon: > python oppgave_f.py 0.75 0.05
"""

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
        N = np.abs(m*g + k*h*(1-l0/np.sqrt(x[i]**2 + h**2)))
        f = -mu*N*np.sign(v[i])
    else:
        f = 0

    a = (-k*x[i]*(1-l0/np.sqrt(x[i]**2 + h**2)) + f)/m
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
ep = -ek
plt.plot(x, ek, label='Kinetic energy')
plt.plot(x, ep, label='Potential energy')
plt.plot(x, ek+ep, label='total energy')
plt.ylabel('Energy [J]')
plt.xlabel('x [m]')
plt.legend(ncol=3)
plt.show()
