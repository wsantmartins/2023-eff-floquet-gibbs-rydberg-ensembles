#script for the operator's definition
import qutip as qt 

#defining the single-body operators
sx = qt.sigmax()
sy = qt.sigmay()
sz = qt.sigmaz()
sp = qt.sigmap()
sm = qt.sigmam()

#defining the two-body operators
n1 = qt.tensor(qt.basis(2,0)*qt.basis(2,0).dag(), qt.qeye(2))
n2 = qt.tensor(qt.qeye(2), qt.basis(2,0)*qt.basis(2,0).dag())
sx1 = qt.tensor(sx, qt.qeye(2))
sy1 = qt.tensor(sy, qt.qeye(2))
sz1 = qt.tensor(sz, qt.qeye(2))
sm1 = qt.tensor(sm, qt.qeye(2))
sp1 = qt.tensor(sp, qt.qeye(2))
sx2 = qt.tensor(qt.qeye(2), sx)
sy2 = qt.tensor(qt.qeye(2), sy)
sz2 = qt.tensor(qt.qeye(2), sz)
sm2 = qt.tensor(qt.qeye(2), sm)
sp2 = qt.tensor(qt.qeye(2), sp)
