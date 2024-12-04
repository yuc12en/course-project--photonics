from eigen_value import *


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
    # dt = get_dt(k0, h0, 0.1, 10)
    # print('time delay: {}'.format(dt))

    # x, E, out = solution_patterns(k0, h0, 10, 10, nf, ns)

    # fig, ax = plt.subplots()
    # ax.plot(x, np.abs(E[0]))
    # ax.plot(x, np.abs(E[2]))
    # plt.show()

    wavelength = 1.5
    k0 = 2*np.pi/wavelength
    nf = 3.5
    ns = 1.5
    nc = 1.5

    mu = 1.257e-6
    w = k0*3e8 * 1e-6
    epsilon=8.85e-12 * 1e6 #F/m * 1e6 um

    h0 = 0.2
    bottom_width = 20-h0/2
    top_width = 20-h0/2
    
    length = 1000
    z_interval = 10
    slope = 0.002
    dx = z_interval*slope

    h_criticle = np.pi*2/np.sqrt(nf**2-ns**2)/k0 + z_interval
    coupling_length = length - (h_criticle-h0)/2/slope  # has to be positive = The end of the waveguide has at least 3 mode. Or there is no intermode coupling
    print(coupling_length)
    # solution_structure(nf, ns, nc)
    # plt.vlines(2*np.pi, 0, 1)
    # plt.show()
    

    x, E, out = solution_patterns(k0, h_criticle, 20-2/h_criticle, 20-2/h_criticle, nf, ns, dx)

    A_mode = 0
    B_mode = 2
    A = E[A_mode]
    B = E[B_mode]
    A = A.astype('complex128')
    B = B.astype('complex128')
    # B = np.zeros(len(B)).astype('complex128')


    yA = out['ys'][A_mode]
    yB = out['ys'][B_mode]
    beta_A = out['beta'][A_mode]
    beta_B = out['beta'][B_mode]


    source_kf = np.pi/10
    source_neff = np.sqrt((source_kf/k0)**2 + ns**2)
    source_ys = k0*np.sqrt(source_neff**2-ns**2)

    x, A, B = intermode_coupling(x, A, B, h0=h_criticle, length=coupling_length, slope=slope, z_interval=z_interval)

    An, C_A = normalized_distribution(x, A, w, mu, beta_A, return_type='coeff')
    Bn, C_B = normalized_distribution(x, B, w, mu, beta_B, return_type='coeff')

    print(C_A, C_B)