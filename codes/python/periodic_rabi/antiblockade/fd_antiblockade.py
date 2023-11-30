import qutip as qt
import matplotlib.pyplot as plt 
import numpy as np
import time
import math
from oprts import*
from stt_qtt import*

start_time = time.time()

delta_values = [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5]
fidelity_ss = []
td_ss = []

for delta in delta_values:
    #periodic detuning
    sw_itm = '0.5 + 0.5*np.sign(np.sin({}*t))'
    sw = sw_itm.format(omega_d)

    #time-independent part of the Hamiltonian for time dependent detuning
    H0 = delta*mag + V*snn
    Hv = magx 

    #Hamiltonians for each half-period for delta variable
    Hp = bigomega*Hv + H0
    Hm = H0


    #Hamiltonian for each half period are defined as H1 and H2
    X = -0.5j*Hp
    Y = -0.5j*Hm

    #first term of the expansion
    P1 = tau*(X + Y)
    
    #P is the only term that has to be calculated explicitly in this program
    def P(i, X, Y, tau):
        P_aux = 0.0
        for k in range(0, i + 1):
            P_aux += ((X**(i-k))*(Y**k))/(math.factorial(i-k)*math.factorial(k))
        return (tau**i)*P_aux

    def Q(i, j):
        if j <= 2 or j >= i:
            return 0.0
        else:
            Qsum = 0.0
            for k in range(1, i - j + 2):
                Qsum += Q(k, 1)*Q(i - j,k - 1)
        return Qsum

    def R(i):
        if i == 1:
            return P1 
        else:
            Qtot = 0.0
            for k in range(2, i + 1):
                Qtot += - Q(i,k)/math.factorial(k)      
            return P(i, X, Y, tau) - Qtot 

    #two particle case: Hamiltonian
    Ht = [H0, [Hv, sw]]

    #operators for calculating the expectation value
    e_ops = [magx, magy, magz]

    #jump operators 
    j_ops = [np.sqrt(gamma)*magm]

    #store the states   
    opts = qt.Options(store_states=True)

    #integrator Lindblad equation
    res_lind = qt.mesolve(Ht, rho0, tspan, j_ops, e_ops, options=opts)

    Rt = 0
    rho_FG = 0
    for k in range(0, Nf):
        Rt += R(k+1)

    H_floquet = 1j*Rt/tau 

    rho_FG = -beta*H_floquet
    rho_FG = rho_FG.expm()
    rho_FG = rho_FG/np.trace(rho_FG) 

    states_at_times = res_lind.states

     

    tr_dist = [qt.tracedist(rho_FG, rho_t) for rho_t in states_at_times]
    fidlty = [qt.fidelity(rho_FG, rho_t) for rho_t in states_at_times]

    fidelity_ss.append(fidlty[nsteps - 1])
    td_ss.append(tr_dist[nsteps - 1])

plt.plot(delta_values, fidelity_ss, marker='o')

print("--- %s seconds ---" % (time.time() - start_time))

plt.show()