import numpy as np
import matplotlib.pyplot as plt

rho = 1.293     # kg/m^3
cd = 1.2
a0 = 0.45       # m^2
w = 0           # m/s
m = 80          # kg
F = 400         # N

dt = 0.01
n = int(np.ceil(10/dt))

t = np.linspace(0, 10, n)
x = np.zeros(n)
v = np.zeros(n)
a = np.zeros(n)
# Alle initialbetingelser er allerede 0 så behøver ikke definere disse.

i_max = 0
for i in range(n-1):
    if x[i-1] > 100:
        break

    a[i+1] = (2*F - rho*cd*a0*(v[i]**2))/(2*m)
    v[i+1] = v[i] + a[i] * dt
    x[i+1] = x[i] + v[i+1] * dt
    i_max += 1

t_max = t[i_max]
t = t[:i_max]; x = x[:i_max]; v = v[:i_max]; a = a[:i_max]

plt.subplot(311)
plt.plot(t, x, 'r', label='x(t)')
plt.ylabel('x [m]')
plt.ylim(0)
plt.xlim(0)

plt.subplot(312)
plt.plot(t, v, 'g', label='v(t)')
plt.ylabel('v [ms^-1]')
plt.ylim(0)
plt.xlim(0)

plt.subplot(313)
plt.plot(t, a, 'b', label='a(t)')
plt.ylabel('a [ms^-2]')
plt.xlabel('t [s]')
plt.xlim(0)

plt.tight_layout()
plt.show()