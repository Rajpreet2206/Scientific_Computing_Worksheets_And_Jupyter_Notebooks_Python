
import numpy as np
import matplotlib.pyplot as plt
#import sympy as sp

def solveODE(n, tau):
    t = np.arange(0,n*tau,tau)

    # analytical solution
    hsint = 0.5*np.sin(t)
    hcost = 0.5*np.cos(t)
    expmt = np.exp(-t)
    avx = np.add(0.5,np.multiply(expmt,np.subtract(hsint,hcost)))
    avy = np.add(-0.5,np.multiply(expmt,np.add(hsint,hcost)))
    
    nvx = np.zeros(n)
    nvy = np.zeros(n)
    erx2 = np.zeros(n)
    ery2 = np.zeros(n)

    # Explicit Euler
    for k in range(1,n):
        nvx[k] = (1-tau)*nvx[k-1]+tau*nvy[k-1]+tau
        nvy[k] = (1-tau)*nvy[k-1]-tau*nvx[k-1]
        erx2[k] = (avx[k]-nvx[k])**2
        ery2[k] = (avy[k]-nvy[k])**2

    er = np.sqrt(np.add(erx2,ery2))

    nsol = {'x':nvx, 'y':nvy}
    asol = {'x':avx, 'y':avy}

    return (t, er, nsol, asol)

#print "sympy version:", sp.__version__

n = 320
tau = 0.03125
t1, er1, nsol1, asol1 = solveODE(n, tau)
t2, er2, nsol2, asol2 = solveODE(int(n/2), 2*tau)
                   
plt.figure()
plt.plot(t1,er1,t2,er2)
plt.title('plot of error compared to analytical solution')
plt.legend(['time step = '+str(tau),'time step = '+str(2*tau)])
plt.xlabel('time')
plt.ylabel('error')

plt.figure()
plt.plot(nsol1['x'], nsol1['y'], nsol2['x'], nsol2['y'], asol1['x'], asol1['y'])
plt.title('numerical solution compared to analytical solution')
plt.legend(['time step = '+str(tau),'time step = '+str(2*tau),'analytical'])
plt.xlabel('vx')
plt.ylabel('vy')

plt.show()