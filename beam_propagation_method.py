import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


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


def gauss_source(x0, h, size=100):
    """
    Generate a beam of gauss profile


    Parameters:
    ----------
    x0: float. characteristic length
    h:  float domain width
    size: int. sampled points

    Returns:
    -------
    E: array. intensity at points
    """
    x = np.linspace(-h/2, h/2, size)
    E = np.exp(-(x/x0)**2)
    return E


def propagate(dz, h, k0, n, E, size=100):
    """
    propagate a beam in the medium. With two phase correction terms

    Parameters:
    ----------
    dz: float. propagation distance
    h: float. domain width
    x0: float.  characteristic length of the beam
    k0: float. wavenumber in vaccum
    n: array of length 2.  Refractive index at the boundary of the interval. 
        Set n[0]=n[1] to get homogeneous medium. Or get linear profile.
    size: int. sampled points size. Must be even

    Returns:
    -------
    newE: array of sampled points size reconstructed at length z.
    """
    # fft components
    n_ave = np.mean(n)
    n = np.linspace(n[0],n[1],size)

    a = h/2/np.pi
    kx = np.concatenate((np.arange(size/2), np.arange(-size/2,0)))/a

    real_beta = n_ave**2*k0**2 - kx**2
    real_beta[real_beta<0] = 0
    phase1 = np.exp(1j*kx**2/(n_ave*k0 +np.sqrt(real_beta))*dz)
    phase2 = np.exp(-1j*(n-n_ave)*k0*dz)

    newE = np.fft.fft(E)*phase1
    newE = np.fft.ifft(newE)*phase2

    return newE


if __name__ == '__main__':
    fig, ax = plt.subplots()

    E = gauss_source(x0=3, h=200, size=512)
    for i in range(10):
        z = (i+1)*10
        newE = propagate(dz=z, h=200, k0=2*np.pi/1, n=[1.5,1.5], size=512, E=E)
        ax.plot(np.arange(len(newE)), np.abs(newE), label='z={}um'.format(z), linewidth=1)

    ax.legend()
    plt.show()
    print(E.max())


