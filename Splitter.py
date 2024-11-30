from Coupler_parameters import *
from eigen_value import *


x, E, out = solution_patterns(wavelength, h0, top_width, bottom_width, nf, ns)

L = 1  # displacement 2L
x1 = x
E1 = E[0]
x2 = x-2*L
E2 = E1

wg_bottom_2 = wg_bottom-2*L
wg_top_2 = wg_top-2*L

E2 = E2[(x2>wg_bottom_2) & (x2<wg_top_2)]
x2 = x2[(x2>wg_bottom_2) & (x2<wg_top_2)]
E1 = E1[(x1>wg_bottom_2) & (x1<wg_top_2)]
x1 = x1[(x1>wg_bottom_2) & (x1<wg_top_2)]

# fig, ax = plt.subplots()
# ax.plot(x1, E1, label='x1')
# ax.plot(x2, E2, label='x2')
# ax.legend()
# plt.show()

kappa = w*mu/4 * (nf**2-ns**2) * simpson(E1*E2, x=x1)
print("Kappa: {}".format(kappa))
print("z0 :{}".format(np.pi/4/kappa))