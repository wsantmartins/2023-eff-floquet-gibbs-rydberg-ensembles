from qutip import* 
import matplotlib.pyplot as plt
import numpy as np 

time = np.linspace(0, 3000,10000)
psi0 = tensor(basis(2,0), basis(2,0))

sx = sigmax()
sy = sigmay()
sz = sigmaz()
sp = sigmap()
sm = sigmam()
n1 = tensor(basis(2,0)*basis(2,0).dag(), qeye(2))
n2 = tensor(qeye(2), basis(2,0)*basis(2,0).dag())

sx1 = tensor(sx, qeye(2))
sy1 = tensor(sy, qeye(2))
sz1 = tensor(sz, qeye(2))
sm1 = tensor(sm, qeye(2))
sp1 = tensor(sp, qeye(2))
sx2 = tensor(qeye(2), sx)
sy2 = tensor(qeye(2), sy)
sz2 = tensor(qeye(2), sz)
sm2 = tensor(qeye(2), sm)
sp2 = tensor(qeye(2), sp)

omega = 1
omega_0 = 1
bigomega = 1
gamma1 = 0.05
V = 1

a_ops = [[sx1 + sx2,"{gamma1} / 2.0 * (w / (2 * pi)) * (w > 0.0)".format(gamma1=gamma1)]]
e_ops = [sx1, sx1 + sx2, sz1, sz1 + sz2]

omega_t = "cos(0.01*t)"
expp = "exp(-1j*t)"
expm = "exp(1j*t)"

Hz = sz1 + sz2
Hp = bigomega*(sp1 + sp2)
Hm = bigomega*(sm1 + sm2)
H_0 = V*(n1 + n2)

H_t = [H_0, [Hz, omega_t], [Hm, expm], [Hp, expp]]

res_brme = brmesolve(H_t, psi0, time, a_ops, e_ops)

plt.figure()

plt.plot(time, res_brme.expect[0], label=r'$\sigma_x$')
plt.plot(time, res_brme.expect[1], label=r'$\sigma^1_x + \sigma^2_x$')
plt.plot(time, res_brme.expect[2], label=r'$\sigma_z$')
plt.plot(time, res_brme.expect[3], label=r'$\sigma^1_z + \sigma^2_z$')

plt.legend()

plt.show()