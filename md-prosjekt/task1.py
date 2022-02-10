import numpy as np
import matplotlib.pyplot as plt

def U(r, eps=1, sigma=1):
    return 4*eps*((sigma/r)**12 - (sigma/r)**6)

def main():
    r = np.linspace(0.9, 3, 10001)
    plt.plot(r, U(r, 1, 1), label='U(r)')
    plt.legend()
    plt.grid()
    plt.xlabel('Distance')
    plt.ylabel('Potential')

if __name__ == '__main__':
    main()
    plt.show()