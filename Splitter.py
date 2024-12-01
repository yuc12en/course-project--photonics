from eigen_value import *

wavelength = 1.5
k0 = 2*np.pi/wavelength
nf = 3.5
ns = 1.5
nc = 1.5

mu = 1.257e-6
w = k0*3e8

def get_z0(L, k0, h0, top_width, bottom_width, nf, ns):
    x, E, out = solution_patterns(k0, h0, top_width, bottom_width, nf, ns)
    x1 = x
    E1 = E[0]
    x2 = x-2*L
    E2 = E1

    wg_bottom = -h0/2
    wg_top = h0/2

    wg_bottom_2 = wg_bottom-2*L
    wg_top_2 = wg_top-2*L

    E2 = E2[(x2>wg_bottom_2) & (x2<wg_top_2)]
    x2 = x2[(x2>wg_bottom_2) & (x2<wg_top_2)]
    E1 = E1[(x1>wg_bottom_2) & (x1<wg_top_2)]
    x1 = x1[(x1>wg_bottom_2) & (x1<wg_top_2)]

    # plt.plot(x2,E2)
    # plt.show()

    kappa = w*mu/4 * (nf**2-ns**2) * simpson(E1*E2, x=x1)
    z0 = np.pi*4/kappa
    return kappa, z0

if __name__ == '__main__':
    L = 1.6
    kappa, z0 = get_z0(L, k0, h0=3, top_width=18.5,bottom_width=18.5, nf=3.5, ns=1.5)
    print(kappa, z0)