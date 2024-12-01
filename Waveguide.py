from eigen_value import *

wavelength = 1.5
k0 = 2*np.pi/wavelength
nf = 3.5
ns = 1.5
nc = 1.5

mu = 1.257e-6
w = k0*3e8

h0 = 3

def get_dt(k0, h0, slope, z_interval):
    length = 1000
    iterations = int(length/z_interval)
    z = np.linspace(0, 1000, iterations)
    neffs = np.zeros(iterations)
    for i in range(iterations):
        h_left = h0 + z_interval*i*slope
        h_right = h0 + z_interval*(i+1)*slope
        h = (h_left+h_right)/2
        out = solve_eigen(k0, h, nf, ns, ns)
        neffs[i] = out['neff'][0]
    
    time_delay = simpson(neffs/3e8, x=z)
    return time_delay

if __name__ == '__main__':
    dt = get_dt(k0, h0, 0.1, 10)
    print('time delay: {}'.format(dt))
