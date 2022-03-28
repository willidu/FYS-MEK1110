import numpy as np
import matplotlib.pyplot as plt
import os

def U(r, eps=1, sigma=1):
    return 4*eps*((sigma/r)**12 - (sigma/r)**6)

def main():
    r = np.linspace(0.9, 3, 10001)
    plt.plot(r, U(r))
    plt.xlabel('Distance')
    plt.ylabel('Potential')
    plt.grid()

if __name__ == '__main__':
    main()
    plt.savefig(os.path.join(os.getcwd(), 'figures/1_a_i.pdf'))
    plt.show()