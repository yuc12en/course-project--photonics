from constants import *
from eigen_value import *
import numpy as np
from matplotlib import pyplot as plt

wavelength = 1
k0 = 2*np.pi/wavelength

def get_dispersion(k0s, h0, length, slope1, z_turn, slope2, z_interval):
    ts = np.zeros(len(k0s))
    for i in range(len(k0s)):
        z, neffs, ts[i] = get_dt(k0s[i], h0, length, slope1, z_turn, slope2, z_interval)
    return ts


def get_dt(k0, h0, length, slope1, z_turn, slope2, z_interval):
    z1 = np.arange(0, z_turn, z_interval)
    z2 = np.arange(z_turn, length + z_interval, z_interval)
    z = np.concatenate((z1,z2))

    h = np.zeros(len(z1)+len(z2))
    h[:len(z1)] = h0+z1*slope1
    h[len(z1):] = h0+z_turn*slope1-(z2-z_turn)*slope2

    neffs = np.zeros(len(h)-1)
    for i in range(len(h)-1):
        h_tmp = (h[i]+h[i+1])/2
        out = solve_eigen(k0, h_tmp, nf, ns, ns)
        neffs[i] = out['neff'][0]

    # use trapz to avoid simpson error
    time_delay = np.trapz(neffs/3e8*1e-6, x=z[:-1])   # x: um * 1e-6 m/um * s/m 

    return z, neffs, time_delay

def intermode_coupling(k0, h0, length, slope1, z_turn, slope2, z_interval, dx):
    w = k0*3e8 *1e6 # um-1 * m/s * um/m = /s  # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m


    z1 = np.arange(0, z_turn, z_interval)
    z2 = np.arange(z_turn, length + z_interval, z_interval)
    z = np.concatenate((z1,z2))
    h = np.zeros(len(z1)+len(z2))
    h[:len(z1)] = h0+z1*slope1
    h[len(z1):] = h0+z_turn*slope1-(z2-z_turn)*slope2

    flag = False
    A = np.nan
    B = np.nan

    Change_CA = []
    Change_CB = []
    Kappas = []
    Ints = []

    for i in range(len(z)-1):
        if z[i] < z_turn:
            h_low = h[i]
            h_high = h[i+1]
            n1 = ns
            n2 = nf
        else:
            h_low = h[i+1]
            h_high = h[i]
            n1 = nf
            n2 = ns

        h_tmp = (h_high+h_low)/2
        out = solve_eigen(k0, h_tmp, nf, ns, ns)

        if len(out['b']) <3:
            continue
        else:
            if flag == False:
                x = np.arange(-20,20,dx)
                A = generate_mode_function(out['kf'][0],out['ys'][0], h_tmp, mode=0)(x)
                B = generate_mode_function(out['kf'][2],out['ys'][2], h_tmp, mode=2)(x)
                beta_A = out['beta'][0]
                beta_B = out['beta'][2]

                An, C_A = normalized_distribution(x, A, w, mu, beta_A, return_type='coeff')
                Bn, C_B = normalized_distribution(x, B, w, mu, beta_B, return_type='coeff')
                C_B = 0
                B = Bn*C_B
                flag = True

        if (z_interval*slope1/2 <= dx*3) or (z_interval*slope2/2 <= dx*3):  # At three x points in the calculation
            raise ValueError("Too sparse x. Increase slope, z_interval or reduce dx")

        indx = (x>h_low/2) & (x<h_high/2)

        x_step = x[indx]
        E1_step = An[indx]
        E2_step = Bn[indx]

        kappa = epsilon*w/4* (n2**2-n1**2) *(np.trapezoid(E1_step*np.conjugate(E2_step), x=x_step))*2
        Ints.append(np.trapezoid(E1_step*np.conjugate(E2_step), x=x_step))
        delta_CB = -2j*kappa/(beta_A-beta_B)*C_A*np.sin((beta_A-beta_B)*z_interval/2)
        delta_CA = -2j*kappa/(beta_B-beta_A)*C_B*np.sin((beta_B-beta_A)*z_interval/2)

        C_A = C_A + delta_CA
        C_B = C_B + delta_CB
        A = An*C_A
        B = Bn*C_B

        Change_CA.append(np.abs(C_A))
        Change_CB.append(np.abs(C_B))
        Kappas.append(kappa)

    return ([[An, C_A], [Bn, C_B]], Change_CA, Change_CB, Kappas, Ints)


def source_coupling_2(k0, length, slope, z_interval ,dx):
    w = k0*3e8 *1e6 # um-1 * m/s * um/m = /s  # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m

    if z_interval*slope/2 <= dx*3:
        raise ValueError("Too sparse x. Increase slope, z_interval or reduce dx")

    h_tip = z_interval*slope
    z = np.arange(0, length+z_interval, z_interval)
    h = z*slope
    # print(h)
    # print(h+h_tip)

    out = solve_eigen(k0, h_tip, nf, ns, ns)
    x = np.arange(-20,20,dx)

    fit_k = 2*np.pi/1.5
    source_kf = 13.244
    source_neff = np.sqrt(nf**2- (source_kf/fit_k)**2)
    source_ys = fit_k*np.sqrt(source_neff**2-ns**2)
    source_beta = fit_k*source_neff

    A = generate_mode_function(source_kf, source_ys, 0.02)(x)
    An, C_A = normalized_distribution(x, A, w, mu, source_beta, return_type='coeff')

    B = generate_mode_function(out['kf'][0],out['ys'][0], h_tip, mode=0)(x)
    beta_t = out['beta'][0]
    Bn, C_B = normalized_distribution(x, B, w, mu, beta_t, return_type='coeff')

    t = beta_t/w/mu/2*simpson(Bn*An, x=x)
    # print(t)

    
    for i in range(len(z)-1):
        h_new = (h[i]+h[i+1])/2
        out = solve_eigen(k0, h_new, nf, ns, ns)

        B_new = generate_mode_function(out['kf'][0], out['ys'][0], h_new)(x)
        beta_t = out['beta'][0]
        Bn_new, C_B = normalized_distribution(x, B_new, w, mu, beta_t, return_type='coeff')

        tmp = beta_t/w/mu/2*np.trapezoid(Bn_new*Bn, x=x)
        if tmp>= 1:
            tmp=1
        t *= tmp
        Bn = Bn_new

    return t


def source_coupling_end(k0, h0, slope, z_interval ,dx):
    w = k0*3e8 *1e6 # um-1 * m/s * um/m = /s  # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m

    if z_interval*slope/2 <= dx*3:
        raise ValueError("Too sparse x. Increase slope, z_interval or reduce dx")

    length = h0/slope

    z = np.arange(0, length+z_interval, z_interval)
    h = h0-z*slope
    # print(h)
    # print(h+h_tip)

    out = solve_eigen(k0, h0, nf, ns, ns)
    x = np.arange(-20,20,dx)
    B = generate_mode_function(out['kf'][0],out['ys'][0], h0, mode=0)(x)
    beta_t = out['beta'][0]
    Bn, C_B = normalized_distribution(x, B, w, mu, beta_t, return_type='coeff')

    t = 1
    for i in range(len(z)-1):
        h_new = (h[i]+h[i+1])/2
        out = solve_eigen(k0, h_new, nf, ns, ns)

        B_new = generate_mode_function(out['kf'][0], out['ys'][0], h_new)(x)
        beta_t = out['beta'][0]
        Bn_new, C_B = normalized_distribution(x, B_new, w, mu, beta_t, return_type='coeff')

        tmp = beta_t/w/mu/2*np.trapezoid(Bn_new*Bn, x=x)
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

    mu = 1.257e-6 *1e-6  #H/m * m/um = H/um
    w = k0*3e8 *1e6 # um-1 * m/s * um/m = /s  # units of mu and w will cancle out. mu: F/m.  w: um-1 * m/s * 10e-6 um/m
    epsilon=8.85e-12 *1e-6 #F/m * m/um  = F/um

    h0 = 0.2
    bottom_width = 20-h0/2
    top_width = 20-h0/2

    # inter-mode coupling(waveguide)

    # slope1 = 0.003
    # slope2 = 0.003
    # z_interval = 10
    # dx = z_interval*np.min([slope1, slope2])/2/4

    # AB, As, Bs, kappas, Ints = intermode_coupling(k0, h0, 1000, slope1, 500, slope2, z_interval, dx)
    # print(Bs[-1]/As[-1])
    # As = [item.max() for item in As]
    # Bs = [item.max() for item in Bs]

    # print(As)
    # print(Bs.shape)
    # print(As[0])
    # print(Bs[0])
    # # Change of the amplititude of higher mode B
    # fig, ax = plt.subplots(2)
    # ax[0].plot(range(len(Bs)), np.abs(Bs), 'r--', label='B')
    # ax[0].set_ylabel('B')
    # ax2 = ax[0].twinx()
    # ax2.plot(range(len(As)), np.abs(As), 'k-', label='A')
    # ax2.set_ylabel('A')
    # ax[0].plot(range(len(As)), np.abs(As), 'k-', label='A')

    # ax[1].plot(range(len(kappas)), np.abs(kappas))
    # ax2 = ax[1].twinx()
    # ax2.scatter(range(len(Ints)), np.abs(Ints))

    # plt.legend()
    # plt.show()

######################################################################################
    # source mode coupling

    # import time
    # start_time = time.time()

    # source-mode courpling(coupler)
    # h0 = 0.2
    # slope = 0.3
    # length = h0/slope
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
    #         ts[i,j] = source_coupling(k0, length[i], slope[j], z_interval, dx)
    
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
    
    ########################3
    # output
    h0 = 0.2
    slope = 3.2 
    slope = 1/3.2
    length = h0/slope
    z_interval = 0.01
    dx = z_interval*slope/2/4
    
    # test_slope = np.arange(0.1,0.40,0.1)
    # ts = np.zeros(len(test_slope))
    # for i in range(len(test_slope)):
    #     dx = z_interval*test_slope[i]/2/4
    #     ts[i] = source_coupling_end(2*np.pi/1.5, 0.2,test_slope[i], z_interval, dx)
    
    # plt.plot(test_slope, ts)
    # plt.show()
    ########################################3
    # plot for transmission
    import time
    start = time.time()
    wavelength_range = np.arange(1.3, 1.7+0.1,0.1)
    t_total = np.zeros(len(wavelength_range))
    for i in range(len(wavelength_range)):
        k0 = np.pi*2/wavelength_range[i]
        t1 = source_coupling_2(k0, length, slope, z_interval, dx)
        t2 = source_coupling_end(k0, 0.2, 0.1, z_interval, z_interval*0.1/2/4)
        t_total[i] = t1*t2
    end = time.time()
    print("Time cost: {}".format(end-start))
    plt.plot(wavelength_range, t_total)
    plt.show()
    
    # plot for time delay









