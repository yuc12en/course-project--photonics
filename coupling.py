from eigen_value import *

import numpy as np
from matplotlib import pyplot as plt

def intermode_coupling(k0, h0, length, slope1, z_turn, slope2, z_interval, dx):
    mu = 1.257e-12
    w = k0*3e8 * 1e-6   # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m

    iterations = int(length/z_interval)

    flag = False
    A = np.nan
    B = np.nan

    Change_CA = []
    Change_CB = []
    for i in range(iterations):
        if z_interval*i <= z_turn:
            h_low = h0 + z_interval*i*slope1
            h_high = h0 + z_interval*(i+1)*slope1
            n1 = ns
            n2 = nf

        if z_interval*i > z_turn:
            h_high = h0 + z_turn*slope1 - (z_interval*i - z_turn)*slope2
            h_low = h0 + z_turn*slope1 - (z_interval*(i+1) - z_turn)*slope2
            n1 = nf
            n2 = ns

        out = solve_eigen(k0, (h_high+h_low)/2, nf, ns, nc)
        h = (h_high+h_low)/2

        if len(out['b']) <3:
            continue
        else:
            if flag == False:
                x = np.arange(-20,20,dx)
                A = generate_mode_function(out['kf'][0],out['ys'][0], h, mode=0)(x)
                B = generate_mode_function(out['kf'][2],out['ys'][2], h, mode=2)(x)
                beta_A = out['beta'][0]
                beta_B = out['beta'][2]
                flag = True

                An, C_A = normalized_distribution(x, A, w, mu, beta_A, return_type='coeff')
                Bn, C_B = normalized_distribution(x, B, w, mu, beta_B, return_type='coeff')
                C_B = 0
                B = Bn*C_B

        if (z_interval*slope1/2 <= dx*3) or (z_interval*slope2/2 <= dx*3):  # At three x points in the calculation
            raise ValueError("Too sparse x. Increase slope, z_interval or reduce dx")

        indx = (x>h_low/2) & (x<h_high/2)

        x_step = x[indx]
        E1_step = An[indx]
        E2_step = Bn[indx]

        kappa = epsilon*w/4* (n1**2-n2**2) *(simpson(E1_step*E2_step, x=x_step))*2

        delta_A = 1j*kappa*C_B*Bn*np.exp(-1j*(beta_A-beta_B)*z_interval/2)*z_interval
        delta_B = 1j* np.conjugate(kappa)*C_A*An *np.exp(1j*(beta_A-beta_B)*z_interval/2)*z_interval
        print(np.abs(delta_B).max())

        # loopnum = 30
        # if i<loopnum:
        #     print(kappa)
            # print(C_A)
            # print(C_B)
            # print('\n')
        # if i==loopnum:
        #     break
        A = C_A*An + delta_A
        B = C_B*Bn + delta_B

        An, C_A = normalized_distribution(x, A, w, mu, beta_A, return_type='coeff')
        Bn, C_B = normalized_distribution(x, B, w, mu, beta_B, return_type='coeff')
        Change_CA.append(np.abs(C_A))
        Change_CB.append(np.abs(C_B))

    return ([[An, C_A], [Bn, C_B]], Change_CA, Change_CB)


def source_coupling(k0, h0, length, slope, z_interval, dx):
    mu = 1.257e-12
    w = k0*3e8 * 1e-6   # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m
    ns = 1.5
    nf = 3.5

    h_tip = z_interval*slope

    iterations = int(length/z_interval)
    if z_interval*slope/2 <= dx*3:
        raise ValueError("Too sparse x. Increase slope, z_interval or reduce dx")

    out = solve_eigen(k0, h_tip, nf, ns, ns)
    x = np.arange(-20,20,dx)

    A = np.exp(-(x/5)**2)
    # beta_A = nf*k0
    # An, C_A = normalized_distribution(x, A, w, mu, beta_A, return_type='coeff')   #  set the beta_i = k0*nf

    B = generate_mode_function(out['kf'][0],out['ys'][0], h_tip, mode=0)(x)
    beta_B = out['beta'][0]
    Bn, C_B = normalized_distribution(x, B, w, mu, beta_B, return_type='coeff')
    C_B = 0
    B = C_B*Bn

    Change_CA = []
    Change_CB = []
    for i in range(iterations):
        h_low = z_interval*i*slope
        h_high = z_interval*(i+1)*slope

        indx = (x>h_low/2) & (x<h_high/2)

        x_step = x[indx]
        E1_step = A[indx]
        E2_step = Bn[indx]


        dBdz = -1j*w/4 *epsilon*(nf**2-ns**2)* simpson(E1_step*E2_step, x=x_step)
        dAdz = -1j*w/4 *epsilon*(nf**2-ns**2)* simpson(E1_step*E1_step, x=x_step)

        A = A + dAdz * z_interval
        B = C_B*Bn + dBdz * z_interval

        loopnum = 3

        Bn, C_B = normalized_distribution(x, B, w, mu, beta_B, return_type='coeff')
        Change_CB.append(C_B)
        Change_CA.append(max(A))

        next_h = z_interval*(2*i+1)*slope/2
        out = solve_eigen(k0, next_h, nf, ns, ns)
        B = generate_mode_function(out['kf'][0],out['ys'][0], next_h, mode=0)(x)
        beta_B = out['beta'][0]
        Bn = normalized_distribution(x, B, w, mu, beta_B)

    return ([[A], [Bn, C_B]], Change_CA, Change_CB)


def source_coupling_2(k0, length, slope, z_interval ,dx):
    mu = 1.257e-12
    w = k0*3e8 * 1e-6   # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m
    ns = 1.5
    nf = 3.5

    iterations = int(length/z_interval)

    if z_interval*slope/2 <= dx*3:
        raise ValueError("Too sparse x. Increase slope, z_interval or reduce dx")

    h_tip = z_interval*slope
    out = solve_eigen(k0, h_tip, nf, ns, ns)

    x = np.arange(-20,20,dx)

    source_kf = np.pi/10
    source_neff = np.sqrt((source_kf/k0)**2 + ns**2)
    source_ys = k0*np.sqrt(source_neff**2-ns**2)
    source_beta = k0*source_neff
    A = generate_mode_function(source_kf, source_ys,5)(x)
    An, C_A = normalized_distribution(x, A, w, mu, source_beta, return_type='coeff')

    B = generate_mode_function(out['kf'][0],out['ys'][0], h_tip, mode=0)(x)
    beta_t = out['beta'][0]
    Bn, C_B = normalized_distribution(x, B, w, mu, beta_t, return_type='coeff')

    t = beta_t/w/mu/2*simpson(Bn*An, x=x)
    

    for i in range(iterations):
        h_low = z_interval*i*slope
        h_high = z_interval*(i+1)*slope
        h_new = (h_low+h_high)/2
        out = solve_eigen(k0, h_new, nf, ns, ns)

        B_new = generate_mode_function(out['kf'][0], out['ys'][0], h_new)(x)
        beta_t = out['beta'][0]
        Bn_new, C_B = normalized_distribution(x, B_new, w, mu, beta_t, return_type='coeff')

        tmp = beta_t/w/mu/2*simpson(Bn_new*Bn, x=x)
        if tmp>= 1:
            tmp=1
        t *= tmp
        Bn = Bn_new

    return t

    


def mode_integration(kf, ys, h, mode='even'):
    if mode=='even' or mode%2==0:
        return 1/ys + 1/np.cos(kf*h/2)**2/(2*kf)*(np.sin(kf*h)+kf*h)
    else:
        return 1/ys + 1/np.sin(kf*h/2)**2/(2*kf)*(kf*h-np.sin(kf*h))



if __name__ == '__main__':
    wavelength = 1.5
    k0 = 2*np.pi/wavelength
    nf = 3.5
    ns = 1.5
    nc = 1.5

    mu = 1.257e-6
    w = k0*3e8
    epsilon=8.85e-12 * 1e6 #F/m * 1e6 um

    h0 = 0.2
    bottom_width = 20-h0/2
    top_width = 20-h0/2

    # inter-mode coupling(waveguide)

    slope1 = 0.003
    slope2 = 0.003
    z_interval = 10
    dx = z_interval*np.min([slope1, slope2])/2/4

    AB, As, Bs = intermode_coupling(k0, h0, 1000, slope1, 500, slope2, z_interval, dx)
    # Change of the amplititude of higher mode B
    fig, ax = plt.subplots()
    ax.plot(range(len(Bs)), np.abs(Bs), 'k--', label='B')
    ax.set_ylabel('B')
    ax2 = ax.twinx()
    ax2.plot(range(len(As)), np.abs(As), 'k-', label='A')
    ax2.set_ylabel('A')

    plt.legend()
    plt.show()

######################################################################################
    # source mode coupling


    # import time
    # start_time = time.time()
    # source-mode courpling(coupler)
    # length = 0.4
    # slope = 0.5
    # z_interval = 0.01
    # dx = z_interval*slope/2/4
    # t = source_coupling_2(k0, length, slope, z_interval, dx)
    # print(t)

    # length = np.linspace(0.2, 1, 10)
    # slope = np.linspace(0.1, 2, 10)
    # z_interval = 0.01
    # ts = np.zeros((len(length), len(slope)))
    # for i in range(len(length)):
    #     for j in range(len(slope)):
    #         dx = z_interval*slope[j]/2/4
    #         ts[i,j] = source_coupling_2(k0, length[i], slope[j], z_interval, dx)
    
    # end_time = time.time()
    # print("Time: {} sec".format(end_time-start_time))

    # fig, ax = plt.subplots()
    # im = ax.imshow(ts)
    # ax.set_xticks(range(len(slope)), labels=np.round(1/slope, 1))
    # ax.set_yticks(range(len(length)), labels=np.round(length, 1))
    # ax.set_xlabel('length/thickness')
    # ax.set_ylabel('length(um)')
    # fig.colorbar(im)
    # plt.show()
    
    # print((ts).max())





