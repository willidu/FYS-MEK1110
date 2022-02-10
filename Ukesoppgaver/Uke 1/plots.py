import matplotlib.pyplot as plt
import numpy as np

A = 3  # m/s^3

def v(t):
    return 1/2 * A * t**2

def a(t):
    return A * t

def x(t):
    return 1/6 * A * t**3 + v(0)

t = np.linspace(0, 5, 101)

plt.plot(t, x(t), label='x(t)')
plt.plot(t, v(t), label='v(t)')
plt.plot(t, a(t), label='a(t)')

plt.xlabel('Time [s]')
plt.legend()
plt.show()