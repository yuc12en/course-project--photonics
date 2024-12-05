from eigen_value import *
from coupling import *

wavelength = 1
k0 = 2*np.pi/wavelength

mu = 1.257e-6 *1e-6  #H/m * m/um = H/um
w = k0*3e8 *1e6 # um-1 * m/s * um/m = /s  # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m
epsilon=8.85e-12 *1e-6 #F/m * m/um  = F/um
nf = 3.5
ns = 1.5

h0 = 0.2
length = 1000

# slope1 = 0.003
# z_turn_frac = 0.5
# z_turn = length*z_turn_frac
# slope2 = slope1 * z_turn_frac/(1-z_turn_frac)
# z_interval = 10
# dx = z_interval*np.min([slope1, slope2])/2/4

# test for dispersion
# wavelength = np.arange(1.300, 1.700,0.02)
# k0s = np.pi*2/wavelength
# ts = get_dispersion(k0s, h0, length, slope1, z_turn, slope2, z_interval)

# ax = solution_structure(nf, ns, ns)
# V1, a = get_normalized_parameters(min(k0s), h0, nf, ns, ns)
# V2, a = get_normalized_parameters(max(k0s), h0, nf, ns, ns)
# ax.vlines([V1,V2], 0, 1, color='k', linestyles='--')
# ax.set_xlim(V1*0.5, V2*2)
# plt.show()

# print(max(ts))
# print(max(ts)-min(ts))

# print(z_turn*slope1 - (length-z_turn)*slope2)

# neff, dt = get_dt(k0, h0, length, slope1, z_turn, slope2, z_interval)
# fig, ax = plt.subplots()
# ax.plot(range(len(neff)), neff)
# plt.show()

######################33
# slope=0.6 testing
# slope1 = 0.002
# z_turn_frac = 0.60000000000000000000001
# z_turn = length*z_turn_frac
# slope2 = slope1 * z_turn_frac/(1-z_turn_frac)
# z_interval = 10
# dx = z_interval*np.min([slope1, slope2])/2/4
# wavelength = np.arange(1.300, 1.700, 0.10)
# k0s = np.pi*2/wavelength
# # dts =  []
# for k in range(len(k0s)):
#     z, neffs, ts = get_dt(k0s[k], h0, length, slope1, z_turn, slope2, z_interval)
    # print(ts)

# neffs = np.load('nerror.npy')
# z = np.load('zerror.npy')

# print(np.trapz(y=neffs/3e8*1e-6, x=z[:-1]))
# plt.plot(z[:-1], neffs)
# plt.show()

########################
# fig, ax = plt.subplots()
# ax.plot(range(len(neffs)), neffs)
# plt.show()

# print(neffs[-1])

# basic solution
# for i in range(iterations):
#     h = h0 - i*h0/iterations
#     sol = solve_eigen(2*np.pi/wavelength, h, nf, ns, ns)
#     print(sol['ys'])

# show the normalized picture and the solution
# ax = solution_structure(nf, ns, ns, n_modes=1)
# plt.show()

# only show normalized picture
# ax = solution_structure(nf, ns, nc)  # n_modes, max_V to adjust the lines and x-axis
# plt.show()

# draw the solution patterns
# fig , ax = plt.subplots()
# for i in range(iterations):
#     h = h0 - i*h0/iterations
#     x, Es, out = solution_patterns(wavelength, h, 3, 3, nf, ns, 0.1, show=False)
#     if i%10 == 0:
#         ax.plot(x, Es[0], label='iterations: {}'.format(i))
# ax.legend()
# plt.show()


# test for BPM
# fig, ax = plt.subplots()
# x, E, out = solution_patterns(wavelength, h0, top_width, bottom_width, nf, ns, dx, show=False)
# E = E[0]
# n = np.ones(len(x))*ns
# n[(x>wg_bottom) & (x<wg_top)]
# ax.plot(x, np.abs(E))
# E = propagate(x, E, nf, k0, 20, 1)
# ax.plot(x, np.abs(E))
# plt.show()