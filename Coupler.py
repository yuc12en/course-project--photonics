from beam_propagation_method import gauss_source, complex_profile, propagate
from eigen_value import solve_eigen, mode_solution

import numpy as np
from matplotlib import pyplot as plt

wavelength = 1.5

# Geometry of waveguide
h0 = 5  
top_width = 3
bottom_width = 3
iterations = 20  # where h=0

# precision
dx = 0.1  # determines the sampled size along x-axis
dz = 10  # um. Stepping size

domain_width = 60  # domain width of BPM. Determing the spacial frequency range after fft range.

# Before coupling
x, Es = mode_solution(wavelength, h0, top_width=top_width, bottom_width=bottom_width, nf=3.5, ns=1.5, dx=dx)
size = x.shape[0]  # sampling size along x-axis
propagated_E = np.zeros((iterations+1, size), dtype='complex128')
propagated_E[0] = Es[0]  # propagate the fundamental mode

# steping in the coupler
for i in range(iterations):
    h = h0 - i*h0/iterations
    neff = solve_eigen(2*np.pi/wavelength, h, nf=3.5, ns=1.5, nc=1.5,return_type='neff')
    neff = neff[0]  # choose fundamental mode
    propagated_E[i+1] = propagate(dz=dz, h=domain_width, k0=2*np.pi/wavelength, n=[neff, neff], size=size, E=propagated_E[i])

fig, ax = plt.subplots()
for i in range(iterations):
    ax.plot(x, propagated_E[i])

plt.show()
