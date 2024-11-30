from bpm import *
from eigen_value import *
from Coupler_parameters import *

import numpy as np
from matplotlib import pyplot as plt


slope = np.linspace(5,10, 20)
ts = np.zeros(len(slope))

for i in range(len(slope)):
    x, E, out = get_inverse_taper_pattern(k0, h0, bottom_width, top_width, slope=slope[i] ,n_interval=int(h0*slope[i]*2), n_step=5, dx=0.1, show=False)

    source = np.exp(-(x/5)**2)
    source = normalized_distribution(x, source, w, mu, k0*nf)
    I = simpson(source*np.conjugate(E), x=x)
    ts[i] = out['beta'][0]/2/w/mu*I

#     out['ax'].set_title('slope:{}'.format(slope[i]))
#     out['ax'].plot(x, source, color='r')
# plt.show()

fig, ax = plt.subplots()
ax.plot(slope, np.abs(ts))
plt.show()


#     ax.plot(x, E, color=plt.get_cmap('viridis')(1-i/len(slope)))
#     ax.plot(x, source, color='r')
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(slope, ts)
# plt.show()