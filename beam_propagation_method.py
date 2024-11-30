import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def overlap(e1, e2, dx):
    I = 0
    for i in range(len(e1)):
        I += dx*e1[i]*e2[i]
    return I


def chop(complex_array, tol=1e-10):
    """
    chop the value to reduce round-off noize
    """
    complex_array.real[np.abs(complex_array.real)<tol] = 0
    complex_array.imag[np.abs(complex_array.imag)<tol] = 0
    return complex_array

def complex_profile(complex_array):
    """
    Exam the complex array
    """
    _, ax = plt.subplots(1,2, figsize=(12,6))
    index = np.arange(len(complex_array))
    ax[0].scatter(index, np.abs(complex_array))
    ax[1].scatter(index, np.angle(complex_array))
    return ax



def propagate(x, E, n, k0, z0, dz):
    """
    propagate a beam in the medium. With two phase correction terms

    Parameters:
    ----------
    x: array. x position of wave
    E: array. Intensity at position x
    n: int of array. refractive index profile
    k0: float. wavenumber in vaccum
    H: float. domain width
    dz: float. stepping width in the BPM. NOT THE PROPAGATION WIDTH OR HOMOGENEOUS LENGTH
    z0: float. propagation length

    Returns:
    -------
    newE: array of sampled points size reconstructed at length z.
    """
    size = len(x)
    H = x[-1] - x[0]
    
    if isinstance(n, int):
        n = np.ones(size)*n
    n_ave = np.mean(n)

    a = H/2/np.pi
    kx = np.concatenate((np.arange(size/2), np.arange(-size/2,0)))/a

    beta_square = n_ave**2*k0**2 - kx**2
    beta_square[beta_square<0] = 0
    phase1 = np.exp(1j*dz*kx**2/(n_ave*k0 +np.sqrt(beta_square)))
    phase2 = np.exp(-1j*(n-n_ave)*k0*dz)

    E = np.fft.fftshift(E)
    iterations = int(z0/dz)
    for i in range(iterations):
        E = np.fft.fft(E)*phase1
        E = np.fft.ifft(E)*phase2
    E = np.fft.fftshift(E)

    return E


if __name__ == '__main__':
    # test for propagate
    # size = 512
    # cladwidth = 200
    # x = np.linspace(-0.5, 0.5, size)*cladwidth
    # E = np.exp(-(x/3)**2)
    # plt.plot(x, np.abs(E))
    # E = propagate(x, E, 1.5, 2*np.pi/1, 800, 4)
    # plt.plot(x, np.abs(E))

    # plt.show()
    pass




