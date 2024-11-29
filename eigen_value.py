import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import root_scalar

def normalized_relation(b, V, a, m):
    left = V*np.sqrt(1-b)
    right = m*np.pi + np.arctan(np.sqrt(b/(1-b))) + np.arctan(np.sqrt((b+a)/(1-b)))
    return left-right

def nomalized_eigen_parameters(k0, h, nf, ns, nc):
    denomenator = np.power(nf, 2)-np.power(ns,2)
    V = k0*h*np.sqrt(denomenator)
    a = (np.power(ns,2)-np.power(nc,2))/denomenator
    return V, a

def solve_eigen(k0, h, nf, ns, nc):

    denomenator = np.power(nf, 2)-np.power(ns,2)
    V = k0*h*np.sqrt(denomenator)
    a = (np.power(ns,2)-np.power(nc,2))/denomenator

    n_allowed_modes = int(np.ceil((V-np.arctan(np.sqrt(a)))/np.pi))
    b = np.zeros(n_allowed_modes)
    for i in range(n_allowed_modes):
        b[i] = root_scalar(normalized_relation, bracket=[0,1-1e-6], args=(V, a, i)).root

    neff = b*denomenator + ns**2
    neff = np.sqrt(neff)
    kf = k0 * np.sqrt(nf**2-neff**2)
    gamma_s = k0 * np.sqrt(neff**2-ns**2)
    gamma_c = k0 * np.sqrt(neff**2-nc**2)
    out = {'b':b}
    out['neff'] = neff
    out['kf'] = kf
    out['gamma_s'] = gamma_s
    out['gamma_c'] = gamma_c

    return out

def overlap(e1, e2, dx):
    I = 0
    for i in range(len(e1)):
        I += dx*e1[i]*e2[i]
    return I


def normalized_distribution(E, w, mu, beta, dx):
    I = overlap(E, E, dx)
    E /= np.sqrt(I)
    E *= np.sqrt(2*w*mu/beta)
    return E

def symmetric_profile(kf, gamma, n_mode, h, bottom_cladding_width, top_cladding_width, dx):
    top = top_cladding_width+h/2
    wg_top = h/2
    wg_bottom = -h/2
    bottom = -bottom_cladding_width-h/2
    x1 = np.arange(wg_top, top, dx)
    x2 = np.arange(wg_bottom, wg_top, dx)
    x3 = np.arange(bottom, wg_bottom, dx)

    E1 = np.exp(-gamma*(x1-h/2))
    if n_mode%2==0:
        E2 = np.cos(kf*x2)/np.cos(kf*h/2)
        E3 = np.exp(gamma*(x3+h/2))
    else:
        E2 = np.sin(kf*x2)/np.sin(kf*h/2)
        E3 = -np.exp(gamma*(x3+h/2))

    E = np.concatenate((E3, E2, E1))
    x = np.concatenate((x3, x2, x1))

    return x, E


def mode_structure(nf, ns, nc, n_modes=5, max_V=20):
    denomenator = np.power(nf, 2)-np.power(ns,2)
    a = (np.power(ns,2)-np.power(nc,2))/denomenator

    fig, ax = plt.subplots()
    for i in range(n_modes):
        start_v = i*np.pi+np.arctan(np.sqrt(a))
        V = np.linspace(start_v, max_V, 100)
        b = np.zeros(len(V))
        for j in range(len(V)):
            b[j] = root_scalar(normalized_relation, bracket=[0,1-1e-10], args=(V[j], a, i)).root
        ax.plot(V, b, label='v={}'.format(i))
    ax.legend()
    ax.set_ylim([0,1])
    return ax


def mode_solution(wavelength, h, top_width=3, bottom_width=3, nf=1.5, ns=1.45, dx=0.1,show=False):
    out = solve_eigen(2*np.pi/wavelength, h, nf, ns, ns)
    kf = out['kf']
    ys = out['gamma_s']

    x = np.arange(-bottom_width-h/2, top_width+h/2, dx)
    xs = np.zeros((len(kf), x.shape[0]))
    Es = np.zeros(xs.shape)

    for i in range(len(kf)):
        x, E = symmetric_profile(kf[i], ys[i], i, h, bottom_cladding_width=top_width, top_cladding_width=bottom_width, dx=dx)
        Es[i] = E

    if show==True:
        fig, ax = plt.subplots(1,len(kf))
        for i in range(len(kf)):
            ax[i].plot(x, Es[i])
        plt.show()

    return x, Es


if __name__ == '__main__':
    nf = 1.5
    ns = 1.45
    nc = 1.45
    h = 5
    wavelength = 1

    out = solve_eigen(2*np.pi/wavelength, h, nf, ns, ns)

    # draw solution
    # b = solve_eigen(2*np.pi/wavelength, h, nf, ns, nc)
    # V, a = normalized_parameter(2*np.pi/wavelength, h, nf, ns, nc)
    # ax = mode_structure(nf, ns, nc)
    # ax.vlines(V, 0, 1, colors='k', linewidth=2, linestyles='-.')
    # ax.scatter(np.ones(len(b))*V, b, facecolor='none', edgecolors='k')

    # draw modes
    # x = np.arange(-5.5, 5.5, 0.1)
    # fig, ax = plt.subplots(1,4)
    # for i in range(len(kf)):
    #     E = symmetric_profile(kf[i], ys[i], i, h, bottom_cladding_width=3, top_cladding_width=3, dx=0.1)
    #     ax[i].plot(x, E)
    # plt.show()

    # draw modes and solution
    # xs, Es = mode_solution(wavelength, h, 3, 3, 1.5, 1.45, 0.1, show=True)


    # ## Coupler
    # wavelength = 1.5
    # nf = 3.5
    # ns = 1.5
    # nc = ns
    # ## Geometry of waveguide
    # h0 = 3
    # top_width = 3
    # bottom_width = 3
    # iterations = 20  # where h=0
    # ## precision
    # dx = 0.1  # determines the sampled size along x-axis
    # dz = 10  # um. Stepping size
    # domain_width = 60  # domain width of BPM. Determing the spacial frequency range after fft range.

    # i=0
    # h = 2
    # # xs, Es = mode_solution(wavelength, h, top_width, bottom_width, 3.5, 1.5, dx, show=True  )
            
    # b = solve_eigen(2*np.pi/wavelength, h, nf, ns, nc)
    # V, a = normalized_parameter(2*np.pi/wavelength, h, nf, ns, nc)
    # ax = mode_structure(nf, ns, nc, n_modes=15, max_V=50)
    # ax.vlines(V, 0, 1, colors='k', linewidth=2, linestyles='-.')
    # ax.scatter(np.ones(len(b))*V, b, facecolor='none', edgecolors='k')
    # plt.show()
