import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import root_scalar

def normalized_relation(b, V, a, m):
    left = V*np.sqrt(1-b)
    right = m*np.pi + np.arctan(np.sqrt(b/(1-b))) + np.arctan(np.sqrt((b+a)/(1-b)))
    return left-right

def solve_eigen(k0, h, nf, ns, nc, mode='wavenumber'):
    if mode == 'wavenumber':
        k0 = k0
    else:
        k0 = 2*np.pi/k0

    denomenator = np.power(nf, 2)-np.power(ns,2)
    V = k0 * h * np.sqrt(denomenator)
    a = (np.power(ns,2)-np.power(nc,2))/denomenator

    n_allowed_modes = np.floor((V-np.arctan(a))/np.pi)
    for i in range(n_allowed_modes):
        b = root_scalar(normalized_relation, args=(V, a, i)).root
    return b

def mode_structure(nf, ns, nc, n_modes =5, mode='wavenumber'):
    denomenator = np.power(nf, 2)-np.power(ns,2)
    a = (np.power(ns,2)-np.power(nc,2))/denomenator

    fig, ax = plt.subplots()
    for i in range(n_modes):
        start_v = i*np.pi+np.arctan(a)
        V = np.linspace(start_v, 20, 100)
        b = np.zeros(len(V))
        for j in range(len(V)):
            b[j] = root_scalar(normalized_relation, x0=0, args=(V[j], a, i)).root
        ax.plot(V, b, label='v={}'.format(i))
    ax.legend()
    return ax


if __name__ == '__main__':
    nf = 1.5
    ns = 1.45
    nc = 1.40
    h = 5
    wavelength = 1
    ax = mode_structure(nf, ns, nc, mode='wavelength')
    plt.show()
    plt.savefig('tmp.png')