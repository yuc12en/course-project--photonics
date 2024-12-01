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

mu = 1.257e-6
w = k0*3e8


# heat map for varing h0 and slope. Find that the bigger h0, the larger t.
# Set 1/y < h0/2 + bottom width to restrict the value of h0
# The optimized slope is around 3.8

slope = np.linspace(2, 10, 10)
h0 = np.linspace(1,6,10)
top_width = 20-h0/2
bottom_width = 3-h0/2
ts = np.zeros((len(h0),len(slope)), dtype='complex')


for i in range(len(h0)):
    for j in range(len(slope)):
        # x, E, out = get_inverse_taper_pattern(k0, h0, bottom_width, top_width, slope=slope[i] ,n_interval=int(h0*slope[i]*2), n_step=5, dx=0.1, show='final')
        x, E, out = get_inverse_taper_pattern(k0, h0[i], bottom_width[i], top_width[i], slope=slope[j] ,n_interval=int(h0[0]*slope[j]*2), n_step=5, dx=0.1, show='none')
        source = np.exp(-(x/5)**2)
        source = normalized_distribution(x, source, w, mu, k0*nf)
        I = simpson(source*np.abs(E), x=x)
        ts[i,j] = out['beta'][0]/2/w/mu * I

fig, ax = plt.subplots()
im = ax.imshow(np.abs(ts))
ax.set_xticks(range(len(slope)), labels=np.round(slope, 1))
ax.set_yticks(range(len(h0)), labels=np.round(h0, 1))
ax.set_xlabel('length/width')
ax.set_ylabel('width(um)')
fig.colorbar(im)
plt.show()


############################################################################################

# It can be shown that the effect of increase h0 is to reduce fft error
# FFT have artifacts Page 271, last paragraph
# one way is to increase the domain size. Which must consider the perturbation behind cladding.
# The other is to apodize the domain. Add attenuation term near the edge
# Neglect the error, it seems the coumpling efficiency is indenedpent from h0

############################################################################################

# Exaime the influence of h0 to coupling efficiency
slope = 3.8
fig, ax = plt.subplots()
for i in range(len(h0)):
    # x, E, out = get_inverse_taper_pattern(k0, h0, bottom_width, top_width, slope=slope[i] ,n_interval=int(h0*slope[i]*2), n_step=5, dx=0.1, show='final')
    x, E, out = get_inverse_taper_pattern(k0, h0[i], bottom_width[i], top_width[i], slope=slope ,n_interval=8, n_step=5, dx=0.1, show='process')
    source = np.exp(-(x/5)**2)
    source = normalized_distribution(x, source, w, mu, k0*nf)
    I = simpson(source*np.abs(E), x=x)
    ts[i] = out['beta'][0]/2/w/mu * I

    # ax.plot(x, np.abs(E), color=plt.get_cmap('viridis')(1-i/len(h0)))
plt.show()


############################################################################################


# Compare with 2 h_tip profile with different h0
# slope = 3.8
# h0 = np.linspace(1,6,5)
# top_width = 20-h0/2
# bottom_width = 3-h0/2
# ts = np.zeros(len(h0))

# x1, E2, out = get_inverse_taper_pattern(k0, h0[0], bottom_width[0], top_width[0], slope=slope ,n_interval=int(h0[0]*slope*2), n_step=5, dx=0.1, show='process')
# source = np.exp(-(x1/5)**2)
# source = normalized_distribution(x1, source, w, mu, k0*nf)
# out['ax'].plot(x1, source)

# x2, E2, out = get_inverse_taper_pattern(k0, h0[-1], bottom_width[-1], top_width[-1], slope=slope ,n_interval=int(h0[0]*slope*2), n_step=5, dx=0.1, show='process')
# source = np.exp(-(x2/5)**2)
# source = normalized_distribution(x2, source, w, mu, k0*nf)
# out['ax'].plot(x2, source)
# plt.show()

############################################################################################

# find that 1/y is far less than h0 or 3. 
# there is no other restriction on h0

# ys = np.zeros(len(h0))
# for i in range(len(h0)):
#     out = solve_eigen(k0, h0[i], nf, ns, ns)
#     ys[i] = out['ys'][0]

# fig, ax = plt.subplots()
# ax.plot(range(len(h0)), 1/ys)
# ax.plot(range(len(h0)), h0/2)
# plt.show()

############################################################################################