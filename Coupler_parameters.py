import numpy as np
# Input 
wavelength = 1.5
k0 = 2*np.pi/wavelength
nf = 3.5
ns = 1.5
nc = 1.5

# Geometry of waveguide
h0 = 5
top_width = 9
bottom_width = 3
wg_top = h0/2
wg_bottom = -h0/2

# Geometry of the taper
slope = 100
dz = 2  # Stepping size
taper_length = h0*slope
iterations = int(taper_length/dz)  # where h=0
h_tip = h0/iterations

# Precision parameters
dx = 0.01  # determines the sampled size along x-axis
size = int((bottom_width+top_width+h0)/dx)
domain_width = 120  # domain width of BPM. Determing the spacial frequency range after fft range.

# Normlization parameters
mu = 1.257e-6
w = k0*3e8
