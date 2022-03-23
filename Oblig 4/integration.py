import numpy as np
import matplotlib.pyplot as plt

u0 = 150
m = 23
x0 = 2
alpha = 39.48

n = 100001
t = np.linspace(0, 6, n)
x = np.zeros_like(t)
dt = t[1] - t[0]

for v_, x_ in ((8.0, -5), (10.0, -5)):
    v = v_      # v0
    x[0] = x_   # x0

    for i in range(n-1):
        a = (-np.sign(x[i])*u0/x0 - alpha*v)/m * (np.abs(x[i])<x0)
        v += a*dt
        x[i+1] = x[i] + v*dt

    plt.plot(t, x, label='x(t)')
    plt.xlabel('Time')
    plt.ylabel('Distance')
    plt.legend()
    plt.savefig(f'plot_integrated_v0_{v_:g}.pdf')
    plt.show()
