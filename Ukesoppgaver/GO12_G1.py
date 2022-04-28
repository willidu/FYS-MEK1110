import numpy as np

R = 3.0  # m
m = 0.50  # kg
v0 = 6.0  # m/s
f = 12/60  # omdreininger / sek
dt = 0.001

t = 0

print('---------------------------------------')
for r0 in ([0, 0, 0], [0, -R, 0]):
    
    print(f'Initial position: {r0}')
    r, v, w = np.asarray([r0, [0, v0, 0], [0, 0, 2*np.pi*f]])  # Initial conditions

    while np.linalg.norm(r) <= R:
        a = -2*np.cross(w, v) - w*np.cross(w, r)
        v += a*dt
        r += v*dt
        t += dt

    x, y, z = r
    print(f'Time: {t:.2f} s')
    print(f'Position (x, y, z): ({x:.2f}, {y:.2f}, {z:.2f})')
    print('---------------------------------------')
