from beam_propagation_method import *
from eigen_value import *

import numpy as np
from matplotlib import pyplot as plt

wavelength = 1.5
k0 = 2*np.pi/wavelength
print('k0: {}'.format(k0))
nf = 3.5
ns = 1.5
nc = 1.5

# Geometry of waveguide
h0 = 10  
top_width = 9
bottom_width = 3
iterations = 60  # where h=0

# precision
dx = 0.01  # determines the sampled size along x-axis
dz = 10  # um. Stepping size

domain_width = 60  # domain width of BPM. Determing the spacial frequency range after fft range.

# Examing the h_tip solution
# print(out)
# V, a =  nomalized_eigen_parameters(k0, h_tip, nf, ns, nc)
# ax = mode_structure(nf, ns, nc, n_modes=5, max_V=20)
# ax.vlines(V, 0, 1, colors='k', linewidth=2, linestyles='-.')
# ax.scatter(np.ones(len(out['b']))*V, out['b'], facecolor='none', edgecolors='k')
# plt.title('Solution')
# plt.show()

n_it = iterations
mu = 1.257e-6
w = k0*3e8

t = np.zeros(n_it)
for i in range(n_it):
    iterations = i+1
    h_tip = h0/iterations
    out = solve_eigen(k0, h_tip, nf, ns, nc)
    x, E = symmetric_profile(out['kf'][0], out['gamma_s'][0], 0, h_tip, bottom_width, top_width, dx)
    E = normalized_distribution(E, w, mu, out['neff'][0]*k0, x[1]-x[0])

    print(2*w*mu/out['neff'][0]/k0)
    print(overlap(E, E, x[1]-x[0]))

    source = np.exp(-x**2/25)
    t[i] = out['neff'][0]*k0 / (k0*3e8 * 1.257e-6) * overlap(E, source, x[1]-x[0])

plt.plot(range(len(t)), t)
plt.show()

plt.plot(x, E)
plt.plot(x, source)
plt.show()



# # Before coupling
# propagated_E = np.zeros((iterations+1, len(E)))
# propagated_E[0] = E  # propagate the fundamental mode

# # steping in the coupler
# for i in range(iterations):
#     h = h0 - i*h0/iterations
#     neff = solve_eigen(2*np.pi/wavelength, h, nf=3.5, ns=1.5, nc=1.5,return_type='neff')
#     neff = neff[0]  # choose fundamental mode
#     propagated_E[i+1] = propagate(dz=dz, h=domain_width, k0=2*np.pi/wavelength, n=[neff, neff], size=size, E=propagated_E[i])

# fig, ax = plt.subplots()
# for i in range(iterations):
#     ax.plot(x, propagated_E[i])

# plt.show()
