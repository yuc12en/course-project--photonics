import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import root_scalar
from scipy.integrate import simpson

def get_normalized_parameters(k0, h, nf, ns, nc):
    denomenator = np.power(nf, 2)-np.power(ns,2)
    V = k0*h*np.sqrt(denomenator)
    a = (np.power(ns,2)-np.power(nc,2))/denomenator
    return V, a

def eigen_problem(b, V, a, m):
    left = V*np.sqrt(1-b)
    right = m*np.pi + np.arctan(np.sqrt(b/(1-b))) + np.arctan(np.sqrt((b+a)/(1-b)))
    return left-right

def solve_eigen(k0, h, nf, ns, nc):

    denomenator = np.power(nf, 2)-np.power(ns,2)
    V = k0*h*np.sqrt(denomenator)
    a = (np.power(ns,2)-np.power(nc,2))/denomenator

    n_allowed_modes = int(np.ceil((V-np.arctan(np.sqrt(a)))/np.pi))
    b = np.zeros(n_allowed_modes)
    for i in range(n_allowed_modes):
        b[i] = root_scalar(eigen_problem, bracket=[0,1-1e-6], args=(V, a, i)).root

    neff = b*denomenator + ns**2
    neff = np.sqrt(neff)
    beta = k0*neff
    kf = k0 * np.sqrt(nf**2-neff**2)
    gamma_s = k0 * np.sqrt(neff**2-ns**2)
    gamma_c = k0 * np.sqrt(neff**2-nc**2)
    out = {'b':b}
    out['beta'] = beta
    out['neff'] = neff
    out['kf'] = kf
    out['ys'] = gamma_s
    out['yc'] = gamma_c
    return out

def get_patterns(out, mode_n, h, bottom_cladding_width, top_cladding_width, dx=0.1):
    top = top_cladding_width+h/2
    bottom = -bottom_cladding_width-h/2
    wg_top = h/2
    wg_bottom = -h/2

    x = np.arange(bottom, top, dx)
    x1 = x[x<wg_bottom]
    x2 = x[(x>=wg_bottom) & (x<=wg_top)]
    x3 = x[x>wg_top]

    y = out['ys'][mode_n]
    kf = out['kf'][mode_n]

    if mode_n%2==0:
        E1 = np.exp(y*(x1+h/2))
        E2 = np.cos(kf*x2)/np.cos(kf*h/2)
    else:
        E1 = -np.exp(y*(x1+h/2))
        E2 = np.sin(kf*x2)/np.sin(kf*h/2)

    E3 = np.exp(-y*(x3-h/2))
    E = np.concatenate((E1, E2, E3))

    return x, E


def solution_structure(nf, ns, nc, n_modes=5, max_V=20, k0=None, h=None):
    denomenator = np.power(nf, 2)-np.power(ns,2)
    a = (np.power(ns,2)-np.power(nc,2))/denomenator

    fig, ax = plt.subplots()

    if k0!= None and h!=None:
        out = solve_eigen(k0, h, nf, ns, nc)
        b = out['b']
        V = k0*h*np.sqrt(denomenator)
        ax.vlines(V, 0, 1, colors='k', linewidth=2, linestyles='--')
        ax.scatter(np.ones(len(b))*V, b,  edgecolors='k')
        max_V = V*1.2
        n_modes = len(b)

    for i in range(n_modes):
        start_v = i*np.pi+np.arctan(np.sqrt(a))
        V = np.linspace(start_v, max_V, 100)
        b = np.zeros(len(V))
        for j in range(len(V)):
            b[j] = root_scalar(eigen_problem, bracket=[0,1-1e-10], args=(V[j], a, i)).root
        ax.plot(V, b, label='v={}'.format(i))
    ax.legend()
    ax.set_ylim([0,1])

    return ax

def normalized_distribution(x, E, w, mu, beta):
    E_square = E*E
    I = simpson(E_square, x=x)
    E /= np.sqrt(I)
    E *= np.sqrt(2*w*mu/beta)
    return E

def solution_patterns(wavelength, h, top_width=3, bottom_width=3, nf=1.5, ns=1.45, dx=0.1,show=False):
    out = solve_eigen(2*np.pi/wavelength, h, nf, ns, ns)
    mode_n = len(out['b'])

    top = top_width+h/2
    bottom = -bottom_width-h/2
    x = np.arange(bottom, top, dx)
    Es = np.zeros((mode_n, len(x)))

    w = 2*np.pi/wavelength * 3e8
    mu = 1.257e-6
    for i in range(mode_n):
        x, E = get_patterns(out, i, h, bottom_cladding_width=bottom_width, top_cladding_width=top_width, dx=dx)
        I = simpson(E, x=x)
        E = normalized_distribution(x, E, w, mu, out['beta'][i])
        Es[i] = E

    if show==True:
        fig, ax = plt.subplots(1,mode_n)
        for i in range(mode_n):
            ax[i].plot(x, Es[i])
        plt.show()

    return x, Es, out


if __name__ == '__main__':

    # basic solution
    # sol = solve_eigen(2*np.pi/wavelength, h, nf, ns, ns)

    # show the normalized picture and the solution
    # ax = solution_structure(nf, ns, nc, k0=np.pi*2/wavelength, h=h)
    # plt.show()

    # only show normalized picture
    # ax = solution_structure(nf, ns, nc)  # n_modes, max_V to adjust the lines and x-axis
    # plt.show()

    # draw the solution patterns
    # ax = solution_patterns(wavelength, h, 3, 3, nf, ns, 0.1, show=True)
    # plt.show()

    pass