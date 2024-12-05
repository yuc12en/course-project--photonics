from eigen_value import *
from coupling import *
import time

start = time.time()

wavelength = 1
k0 = 2*np.pi/wavelength
mu = 1.257e-6 *1e-6  #H/m * m/um = H/um
w = k0*3e8 *1e6 # um-1 * m/s * um/m = /s  # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m
epsilon=8.85e-12 *1e-6 #F/m * m/um  = F/um
nf = 3.5
ns = 1.5


# Coupling B/A
h0 = 0.2
bottom_width = 20-h0/2
top_width = 20-h0/2
length = 1000
slope1 = np.arange(0.001, 0.010, 0.001)
z_turn_frac = np.arange(0.2, 0.8, 0.1)
wavelength = np.arange(1.300, 1.700, 0.10)

# inter-mode coupling(waveguide)
#####################################################################################
z_interval = 10

ts = np.zeros((len(slope1), len(z_turn_frac)))
for i in range(len(slope1)):
    for j in range(len(z_turn_frac)):
        z_turn = length*z_turn_frac[j]
        slope2 = slope1[i] * z_turn_frac[j]/(1-z_turn_frac[j])
        dx = z_interval*np.min([slope1[i], slope2])/2/4
        AB, As, Bs, kappas, Ints = intermode_coupling(k0, h0, length, slope1[i], z_turn, slope2, z_interval, dx)
        ts[i,j] = Bs[-1]/As[-1]

end = time.time()


end = time.time()
print("Total Iterations: {}".format(len(slope1)*len(z_turn_frac)))
print("Total time cost: {}".format(end-start))

fig, ax = plt.subplots()
im = ax.imshow(ts)
ax.set_xticks(range(len(z_turn_frac)), labels=np.round(z_turn_frac,1))
ax.set_yticks(range(len(slope1)), labels=np.round(slope1,3))
ax.set_xlabel('z_turn_frac')
ax.set_ylabel('slope')
fig.colorbar(im)
plt.show()

#############################################################
# dispersion
# slope1 = np.arange(0.001, 0.006, 0.001)
# z_turn_frac = np.arange(0.2, 0.8, 0.2)
# wavelength = np.arange(1.300, 1.700, 0.10)

# k0s = np.pi*2/wavelength
# ts = np.zeros((len(slope1), len(z_turn_frac)))
# for i in range(len(slope1)):
#     for j in range(len(z_turn_frac)):
#         z_turn = length*z_turn_frac[j]
#         slope2 = slope1[i] * z_turn_frac[j]/(1-z_turn_frac[j])
#         z_interval = 10
#         dx = z_interval*np.min([slope1[i], slope2])/2/4
#         t_tmp = get_dispersion(k0s, h0, length, slope1[i], z_turn, slope2, z_interval)
#         ts[i,j] = np.ptp(t_tmp)

# # # neffs, dt = get_dt(k0, h0, length, slope1, z_turn, slope2, z_interval)
# fig, ax = plt.subplots()
# im = ax.imshow(ts)
# ax.set_xticks(range(len(z_turn_frac)), labels=np.round(z_turn_frac,1))
# ax.set_yticks(range(len(slope1)), labels=slope1)
# ax.set_xlabel('z_turn_frac')
# ax.set_ylabel('slope')
# fig.colorbar(im)
# plt.show()
    
##################################################################################
# Examine the effect of z_tunr_frac
# slope1 = 0.002
# z_interval = 10
# fig, ax = plt.subplots(1,4)
# for j in range(len(z_turn_frac)):
#     z_turn = length*z_turn_frac[j]
#     slope2 = slope1 * z_turn_frac[j]/(1-z_turn_frac[j])

#     ax[j].set_title('{}'.format(z_turn_frac[j]))
#     dts = []

#     for k in range(len(k0s)):
#         z, neffs, ts = get_dt(k0s[k], h0, length, slope1, z_turn, slope2, z_interval)
#         dts.append(ts)
#         ax[j].plot(z[:-1], neffs, color=plt.get_cmap('viridis')(1-k/len(k0s)), label=k0s[k])
#     ax[j].text(0.1, 0.5, 'Dispersion: {:.3g}'.format(max(dts)-min(dts)), transform=ax[j].transAxes)
#     t_tmp = get_dispersion(k0s, h0, length, slope1, z_turn, slope2, z_interval)
# plt.show()

# end = time.time()
# print("Total Iterations: {}".format(len(k0s)*len(z_turn_frac)))
# print("Total time cost: {}".format(end-start))

# plt.show()
# fig, ax = plt.subplots()
# im = ax.imshow(ts)
# ax.set_xticks(range(len(z_turn_frac)), labels=z_turn_frac)
# ax.set_yticks(range(len(slope1)), labels=slope1)
# ax.set_xlabel('length/thickness')
# ax.set_ylabel('length(um)')
# fig.colorbar(im)
# plt.show()


