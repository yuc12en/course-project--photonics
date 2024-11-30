from eigen_value import *
from Coupler_parameters import *
from beam_propagation_method import *

nf = 1.5
ns = 1.45
nc = 1.45
h = 5
wavelength = 1


# basic solution
# for i in range(iterations):
#     h = h0 - i*h0/iterations
#     sol = solve_eigen(2*np.pi/wavelength, h, nf, ns, ns)
#     print(sol['ys'])

# show the normalized picture and the solution
# ax = solution_structure(nf, ns, nc, k0=np.pi*2/wavelength, h=h)
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