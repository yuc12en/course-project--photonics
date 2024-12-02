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

def intermode_coupling(x, A, B, h0, slope, z_interval):

    length = 1000
    iterations = int(length/z_interval)

    for i in range(iterations):
        An, C_A = normalized_distribution(x, A, w, mu, beta_A, return_type='coeff')
        Bn, C_B = normalized_distribution(x, B, w, mu, beta_B, return_type='coeff')
        if i==0:
            C_B = 0
            B = Bn*C_B

        h_left = h0 + z_interval*i*slope
        h_right = h0 + z_interval*(i+1)*slope

        indx1 = (x>h_left/2) & (x<h_right)
        indx2 = (x<-h_left/2) & (x>-h_right)

        if sum(indx1) <3:
            return x, A, B

        x_top = x[indx1]
        x_bottom = x[indx2]
        E1_top = An[indx1]
        E1_bottom = An[indx2]
        E2_top = Bn[indx1]
        E2_bottom = Bn[indx2]

        kappa = epsilon*w/4* (ns**2-nf**2) *(simpson(E1_top*E2_top, x=x_top) + simpson(E1_bottom*E2_bottom, x=x_bottom))

        delta_A = 1j*kappa*B*np.exp(-1j*(beta_A-beta_B)*z_interval/2)*z_interval
        delta_B = 1j*np.conjugate(kappa)*A*np.exp(1j*(beta_A-beta_B)*z_interval/2)*z_interval
        # print(delta_A)
        # print(delta_B)
        # print('\n\n')
        # if i==3:
        #     break

        A = C_A*An + delta_A
        B = C_B*Bn + delta_B

        
    return x, A, B

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
    w = k0*3e8
    epsilon=8.85e-12 * 1e6 #F/m * 1e6 um

    h0 = 3
    bottom_width = 9
    top_width = 9

    
    z_interval = 10
    slope = 0.002
    dx = z_interval*slope
    x, E, out = solution_patterns(k0, h0, bottom_width, top_width, nf, ns, dx)

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


    x, A, B = intermode_coupling(x, A, B, h0=h0, slope=slope, z_interval=z_interval)

    An, C_A = normalized_distribution(x, A, w, mu, beta_A, return_type='coeff')
    Bn, C_B = normalized_distribution(x, B, w, mu, beta_B, return_type='coeff')

    print(C_A, C_B)