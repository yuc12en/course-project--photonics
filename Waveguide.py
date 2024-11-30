from eigen_value import *
from Coupler_parameters import *

h0 = 5
slope = 0.1

length = 1000
# dz = 1  # z position precision
z_interval = 10  # discrete precision
# z = np.arange(0,length, dz)

iterations = int(length/z_interval)

z = np.linspace(0,1000, iterations)
neffs = np.zeros(iterations)

for i in range(iterations):
    h_left = h0 + z_interval*i*slope
    h_right = h0 + z_interval*(i+1)*slope
    h = (h_left+h_right)/2
    out = solve_eigen(k0, h, nf, ns, ns)
    neffs[i] = out['neff'][0]


time_delay = simpson(neffs/3e8, x=z)
print("Wavelength: {},\tTime delay:{}".format(wavelength,time_delay))

# fig, ax = plt.subplots()
# ax. plot(z, neffs)
# plt.show()