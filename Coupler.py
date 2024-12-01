from bpm import *
from eigen_value import *

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
# Parameters setting
wavelength = 1.5
k0 = 2*np.pi/wavelength
nf = 3.5
ns = 1.5
nc = 1.5


h0 = 2  # h0 < 6
top_width = 20
bottom_width = 3-h0/2

mu = 1.257e-6
w = k0*3e8


slope = np.linspace(5, 6, 3)
h0 = np.linspace(1,6,10)
ts = np.zeros((len(h0),len(slope)), dtype='complex')

# heat map for varing h0 and slope. Find that the bigger h0, the larger t.
# Set 1/y < h0/2 + bottom width to restrict the value of h0
# The optimized slope is around 6
# for i in range(len(h0)):
#     for j in range(len(slope)):
#         x, E, out = get_inverse_taper_pattern(k0, h0, bottom_width, top_width, slope=slope[i] ,n_interval=int(h0*slope[i]*2), n_step=5, dx=0.1, show='final')
#         # x, E, out = get_inverse_taper_pattern(k0, h0[i], bottom_width, top_width, slope=slope[j] ,n_interval=int(h0[i]*slope[j]*2), n_step=5, dx=0.1, show='final')
#         source = np.exp(-(x/5)**2)
#         source = normalized_distribution(x, source, w, mu, k0*nf)
#         I = simpson(source*np.abs(E), x=x)
#         ts[i,j] = out['beta'][0]/2/w/mu * I

# fig, ax = plt.subplots()
# im = ax.imshow(np.abs(ts))
# ax.set_xticks(range(len(slope)), labels=slope)
# ax.set_yticks(range(len(h0)), labels=np.round(h0, 1))
# ax.set_xlabel('length/width')
# ax.set_ylabel('width(um)')
# fig.colorbar(im)
# plt.show()


# find that 1/y is far less than h0 or 3. Must count for the perturbated distribution of cladding.
# ys = np.zeros(len(h0))
# for i in range(len(h0)):
#     out = solve_eigen(k0, h0[i], nf, ns, ns)
#     ys[i] = out['ys'][0]

# fig, ax = plt.subplots()
# ax.plot(range(len(h0)), 1/ys)
# ax.plot(range(len(h0)), h0/2)
# plt.show()

# Found that the behavior may be borught by the artifact form FFT
# one way is to increase the domain size. Which must consider the perturbation behind cladding.
# The other is to apodize the domain. Add attenuation term near the edge
# slope = 6
# h0 = np.linspace(1,6,5)
# ts = np.zeros(len(h0), dtype='complex')

# for i in range(len(h0)):
#     x, E, out = get_inverse_taper_pattern(k0, h0[i], bottom_width, top_width, slope=slope,n_interval=int(h0[i]*slope*2), n_step=5, dx=0.1, show='process')
#     source = np.exp(-(x/5)**2)
#     source = normalized_distribution(x, source, w, mu, k0*nf)
#     I = simpson(source*np.abs(E), x=x)
#     ts[i] = out['beta'][0]/2/w/mu * I
#     out['ax'].plot(x, source, color='r')
#     out['ax'].set_title('h0: {}um'.format(h0[i]))
# plt.show()
