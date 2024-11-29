from beam_propagation_method import *
from eigen_value import *
from Coupler_parameters import *

import numpy as np
from matplotlib import pyplot as plt


# Fundamental mode of waveguide
x, Es, out = solution_patterns(wavelength, h0, top_width, bottom_width, nf, ns, dx, show=False)
propagated_E = np.zeros((iterations+1, size), dtype='complex')
propagated_E[0] = Es[0]  # fundamental mode

# fig, ax = plt.subplots()
# ax.plot(x, Es[0])
# plt.show()

# steping in the coupler 
# Notice: 
# 1. the domain width changes with propagation. While it negelcts in the simulation
for i in range(iterations):
    h = h0 - i*h0/iterations
    neff = solve_eigen(k0, h, nf, ns, nc)['neff'][0]
    propagated_E[i+1] = propagate(dz=dz, h=domain_width, k0=k0, n=[neff, neff], size=size,E=propagated_E[i])

fig, ax = plt.subplots()
for i in range(iterations):
    if i%int(iterations/3) == 0:
        ax.plot(x, np.abs(propagated_E[i]), c='k')
    else:
        continue

source = np.exp(-x**2/25)
ax.plot(x, source, color='r')
plt.show()
