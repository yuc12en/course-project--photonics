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


def BPM(z, h, x0, k, n,size=100):
    """

    Parameters:
    ----------
    z: float. propagation distance
    h: float. domain width
    x0: characteristic length of the beam
    k: wave vector
    n: refractive index
    size: int. sampled points size

    Returns:
    -------
    newE: array of sampled points size reconstructed at length z.
    """
    dx = h/size
    x = np.arange(-h/2, h/2, dx)
    E = np.exp(-x**2/x0)
    Ef = np.fft.fft(E)
    Ef = chop(Ef)

    wavelength = 2*np.pi/(k*n)
    s = np.arange(len(x)/2)
    s = np.concatenate((s, s[::-1]+1))
    s = s**2
    phase_shift = k*n * z * (1 + 0.5*(wavelength/h)**2* s)

    newEf = Ef*np.exp(1j*phase_shift)
    newEf = np.abs(Ef) * np.exp(1j * (np.angle(Ef)+phase_shift))
    newE = np.fft.ifft(newEf)
    newE = chop(newE)

    return newE


if __name__ == '__main__':
    fig, ax = plt.subplots()

    for i in range(5):
        z = 200 + i*100
        newE = BPM(z, 120, 8, 2*np.pi/0.8, 1)
        ax.scatter(np.arange(len(newE)), np.abs(newE), label='z={}um'.format(z), s=1)

    ax.legend()
    plt.show()


