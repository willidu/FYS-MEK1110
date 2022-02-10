# GF oppgave 3.3

import numpy as np
import matplotlib.pyplot as plt

# a)
x = np.linspace(-2*np.pi, 2*np.pi, 100)
y = lambda x: 4*np.sin(2*x)
plt.plot(x, y(x), label='4sin(2x)')
plt.legend(loc=1)
plt.show()

# b)
x = np.linspace(-10, 2, 100)
y = lambda x: x**2*np.exp(0.5*x)
plt.plot(x, y(x), label='x*2*e^(0.5x')
plt.legend()
plt.show()

# c)
t = np.linspace(0, 8*np.pi, 1000)
x = lambda t: t*np.cos(t)
y = lambda t: t*np.sin(t)
plt.plot(x(t), y(t), label='x=tcos(t), y=tsin(t)')
plt.legend()
plt.show()

# d)
x = np.linspace(-np.pi, np.pi, 100)
y = x.copy()
xx, yy = np.meshgrid(x, y)
z = lambda x, y: np.sin(x)*np.sin(y)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(xx, yy, z(xx, yy), rstride=1, cstride=1, cmap=plt.cm.jet,
linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("h")
plt.show()

# e)
x = np.linspace(1, 5, 401)
y = x.copy()
xx, yy = np.meshgrid(x, y)

#z = np.log(xx*yy) / (xx**2 + yy**2)
z = lambda x, y: np.log(x*y) / (x**2 + y**2)
contour = plt.contour(xx, yy, z(xx, yy))
plt.clabel(contour)
plt.show()

# f) 
t = np.linspace(-5, 5, 11)
x, y = np.meshgrid(t, t)
vx = x / (x**2 + y**2)
vy = y / (x**2 + y**2)

for a, b in zip((vx, vx, vy, vy), (vy, -vy, vx, -vx)):
    plt.quiver(x, y, a, b, scale=20)
    plt.axis('equal')
    plt.show()