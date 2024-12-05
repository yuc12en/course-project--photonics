import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import root_scalar
from scipy.integrate import simpson

def overlap(e1, e2, dx):
    I = 0
    for i in range(len(e1)):
        I += dx*e1[i]*e2[i]
    return I

def get_normalized_parameters(k0, h, nf, ns, nc):
    denomenator = np.power(nf, 2)-np.power(ns,2)
    V = k0*h*np.sqrt(denomenator)
    a = (np.power(ns,2)-np.power(nc,2))/denomenator
    return V, a

def eigen_problem(b, V, a, m):
    left = V*np.sqrt(1-b)
    right = m*np.pi + np.arctan(np.sqrt(b/(1-b))) + np.arctan(np.sqrt((b+a)/(1-b)))
    return left-right


# 1
def solve_eigen(k0, h, nf, ns, nc):

    denomenator = nf**2-ns**2
    V = k0*h*np.sqrt(denomenator)
    a = (ns**2-nc**2)/denomenator

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

def generate_mode_function(kf, ys, h, mode='even'):
    if mode=='even' or mode%2==0:
        def E_n(x):
            x1 = x[x<-h/2]
            x2 = x[(x>=-h/2) & (x<=h/2)]
            x3 = x[x>h/2]
            E1 = np.exp(ys*(x1+h/2))
            E2 = np.cos(kf*x2)/np.cos(kf*h/2)
            E3 = np.exp(-ys*(x3-h/2))
            return np.concatenate((E1,E2,E3))
    else:
        def E_n(x):
            x1 = x[x<-h/2]
            x2 = x[(x>=-h/2) & (x<=h/2)]
            x3 = x[x>h/2]
            E1 = -np.exp(ys*(x1+h/2))
            E2 = np.sin(kf*x2)/np.sin(kf*h/2)
            E3 = np.exp(-ys*(x3-h/2))
            return np.concatenate((E1,E2,E3))

    return E_n
# 4

# 2
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
    # ax.set_ylim([0,1])

    return ax

def normalized_distribution(x, E, w, mu, beta, return_type='none'):
    E_square = E*np.conjugate(E)
    I = simpson(E_square, x=x)
    E /= np.sqrt(I)
    E *= np.sqrt(2*w*mu/beta)
    if return_type == 'coeff':
        return E, np.sqrt(I)/np.sqrt(2*w*mu/beta)
    else:
        return E

# 3
def solution_patterns(k0, h, bottom_width=3, top_width=3, nf=3.5, ns=1.5, dx=0.1,show=False):
    out = solve_eigen(k0, h, nf, ns, ns)
    mode_n = len(out['b'])

    top = top_width+h/2
    bottom = -bottom_width-h/2
    x = np.arange(bottom, top, dx)
    Es = np.ones((mode_n, len(x)))

    mu = 1.257e-6 *1e-6  #H/m * m/um = H/um
    w = k0*3e8 *1e6 # um-1 * m/s * um/m = /s  # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m
    epsilon=8.85e-12 *1e-6 #F/m * m/um  = F/um
    for i in range(mode_n):
        E_fun = generate_mode_function(out['kf'][i], out['ys'][i], h, i)
        E = E_fun(x)
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
    wavelength = 1.5
    k0 = 2*np.pi/wavelength
    nf = 3.5
    ns = 1.5
    nc = 1.5

    mu = 1.257e-6 *1e-6  #H/m * m/um = H/um
    w = k0*3e8 *1e6 # um-1 * m/s * um/m = /s  # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m
    epsilon=8.85e-12 *1e-6 #F/m * m/um  = F/um

    # basic solution
    # sol = solve_eigen(2*np.pi/wavelength, h, nf, ns, ns)

    # show the normalized picture and the solution
    ax = solution_structure(nf, ns, nc, k0=np.pi*2/wavelength, h=1.2)
    plt.show()

    # only show normalized picture
    # ax = solution_structure(nf, ns, nc)  # n_modes, max_V to adjust the lines and x-axis
    # plt.show()

    # draw the solution patterns
    # ax = solution_patterns(wavelength, 2.5, 3, 3, nf, ns, 0.1, show=True)
    # plt.show()

    pass