from eigen_value import *


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

    kappa = w*mu/4 * (nf**2-ns**2) * simpson(E1*E2, x=x1)
    z0 = np.pi/kappa/4
    return kappa, z0

def allowed_bending_radius(k0, nf, ns):
    # max of confined single mode
    V = np.pi
    h = V/np.sqrt(nf**2-ns**2)/k0
    h = np.array([h, h/10])
    r = np.array([5, 20])

    return h, r



if __name__ == '__main__':

    # Kappa and z0
    wavelength = 1.5
    k0 = 2*np.pi/wavelength
    nf = 3.5
    ns = 1.5
    nc = 1.5

    mu = 1.257e-6
    w = k0*3e8
    
    # h0 = 0.2
    # distance = 0.1
    # L0 = h0+0.1
    # kappa, z0 = get_z0(L0, k0, h0, top_width=20-h0/2,bottom_width=20-h0/2, nf=3.5, ns=1.5)
    # print(kappa, z0)

    ##########################################################################################
    # Bending radius
    # h, r = allowed_bending_radius(k0, nf, ns)
    # p = np.polyfit(h, r, deg=1)
    # r0 = np.polyval(p, h0)
    # print(r0)

    # fig, ax = plt.subplots()
    # ax.plot(h,r)
    # ax.scatter(h0, r0)
    # ax.set_xlabel('Thickness (um)')
    # ax.set_ylabel('Minimum bending radius (um)')
    # plt.show()

    # print("Minmum Bending width for h0={:.2f} is {:.2f}".format(h0, np.polyval(p, h0)))

    # plt.show()

    ###########################################
    ## sensitive to z but not sensitive to L ##
    ###########################################

    # sensitivity to L
    # h0 = 0.2
    # L0 = h0 + 0.1
    # fig, ax = plt.subplots()
    # As = np.zeros(20, dtype='complex128')
    # Bs = np.zeros(20, dtype='complex128')
    # kappa, z0 = get_z0(L0, k0, h0, top_width=20-h0/2, bottom_width=20-h0/2, nf=3.5, ns=1.5)

    # for i in range(20):
    #     L = L0 + 0.001*i
    #     kappa, _ = get_z0(L, k0, h0, top_width=20-h0/2, bottom_width=20-h0/2, nf=3.5, ns=1.5)
    #     As[i] = np.cos(kappa*z0)
    #     Bs[i] = -1j*np.sin(kappa*z0)

    # ax.scatter(range(20), As, color='k', s=10, label='A')
    # ax.scatter(range(20), np.abs(Bs), color='k', facecolor='None', s=100, label='B')
    # ax.legend()
    # plt.show()


    # sensitive to z
    h0 = 0.2
    L0 = h0 + 0.3
    fig, ax = plt.subplots()
    As = np.zeros(20, dtype='complex128')
    Bs = np.zeros(20, dtype='complex128')
    kappa, z0 = get_z0(L0, k0, h0, top_width=20-h0/2, bottom_width=20-h0/2, nf=3.5, ns=1.5)
    print(kappa, z0)

    zs = np.linspace(z0, z0+0.005, 100)
    As = np.cos(kappa*zs)
    Bs = -1j*np.sin(kappa*zs)

    ax.scatter(zs-z0, np.abs(As), color='k', s=10, label='A')
    ax.scatter(zs-z0, np.abs(Bs), color='k', facecolor='None', s=100, label='B')
    ax.set_ylim([0,1])
    ax.legend()
    plt.show()


