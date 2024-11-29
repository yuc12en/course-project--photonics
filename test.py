from eigen_value import *

nf = 1.5
ns = 1.45
nc = 1.45
h = 5
wavelength = 1

# basic solution
# sol = solve_eigen(2*np.pi/wavelength, h, nf, ns, ns)

# show the normalized picture and the solution
# ax = solution_structure(nf, ns, nc, k0=np.pi*2/wavelength, h=h)
# plt.show()

# only show normalized picture
# ax = solution_structure(nf, ns, nc)  # n_modes, max_V to adjust the lines and x-axis
# plt.show()

# draw the solution patterns
# x, Es, out = solution_patterns(wavelength, h, 3, 3, nf, ns, 0.1, show=True)
# plt.show()
