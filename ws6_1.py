
import numpy as np
import matplotlib.pyplot as plt
import time

# initial setup
timesteps = 1000; N = 20; tau = 0.001; Dcoeff = 1; h = 1./N
leftBoundary = 0; rightBoundary = 1
T = np.zeros(N+1); rhsvec = np.zeros(N+1)
T[0]=leftBoundary
T[N]=rightBoundary
rhsvec[0]=leftBoundary
rhsvec[N]=rightBoundary
# plot the inital

plt.ion()

x = np.linspace(0.,1.,N+1)

# time loop
coeff = Dcoeff/h**2
for t in range(0,timesteps):

    if t%20 == 0:
    	plt.plot(x, T)
    	plt.draw()
    	time.sleep(0.05)
    
    rhsvec[1:N] = [coeff*(T[k-1] - 2*T[k] + T[k+1]) for k in range(1,N)] # +b is included in the boundary terms (k=1, k=N-1)
    T[1:N] = [T[k] + tau*rhsvec[k] + coeff* (tau**2/2) * (rhsvec[k-1] - 2*rhsvec[k] + rhsvec[k+1]) for k in range(1,N)] # A*rhs is properly calculated (since boundary terms are fixed

plt.show(block=True)
