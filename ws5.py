import numpy as np
import matplotlib.pyplot as plt
from pylab import *

def plotDirectionField(X, Y, U, V, pname):
    """ simple plotting of direction field."""
    M = sqrt(pow(U, 2) + pow(V, 2))
    U = np.divide(U,M)
    V = np.divide(V,M)
    #Q = quiver( X, Y, U, V, M, units='xy', pivot='tail', width=0.010, scale=4.0)
    Q = quiver( X, Y, U, V, units='xy', pivot='tail', width=0.008, scale=5.0)
    title(pname)

############################################
################### Task a) ################
# initial setup
timesteps = 4; tau = 1.
p = np.zeros(timesteps+1)
time = [k*tau for k in range(0,timesteps+1)]

# Verhust - Logistic Growth model
alpha = 1./2; beta = 1./2
#f = lambda t, p: (alpha - beta*p)*p
# Verhust - saturation
f = lambda t, p: alpha - beta*p

p[0] = 0.5
fold = f(0,p[0])
print("f[0]", fold)
p[1] = p[0] + tau*fold

# second order
# ptemp = p[0] + tau*f(0,p[0])/2.
# p[1] = p[0] + tau*f(tau/2., ptemp)

# exact solution value
# p[1] = (p[0] + alpha/beta*(exp(beta*tau)-1))*exp(-beta*tau)


for n in range(1,timesteps):
    fnew = f(n*tau,p[n])
    print("time step index", n)
    print("fnew", fnew)
    p[n+1] = 5*p[n-1] - 4*p[n] + tau*(2*fold + 4*fnew)
    fold = fnew

print("task a) numerical scheme output:", p)

# plotting ranges
dt = .25; dp = .1
t0 = 0; t1 = timesteps*tau+dt/2.
p0 = 0.0; p1 = 3.1

plt.figure()

T,P = meshgrid( arange(-t0,t1+dt,dt),arange(p0,p1+dp,dp) )
U = np.array([ [1.]*len(Pr) for Pr in P ])
V = f(0,P)
plotDirectionField(T, P, U, V, "")
# draw horizontal lines at the timestep marks
for t in range(0,timesteps+1):
    plt.plot([t*tau, t*tau], [p0, p1], 'k--')
plt.plot(time,p,'bo',time,p,linewidth=3)
plt.xlim((t0,t1))
plt.ylim((p0,p1))
plt.xlabel('t')
plt.ylabel('p(t)')


############################################
################### Task b) ################

# initial setup
timesteps = 4; tau = 1.
p = np.zeros(timesteps+1)
time = [k*tau for k in range(0,timesteps+1)]


p[0] = 3.
fold = f(0,p[0])
#print "f[0]", fold
#p[1] = p[0] + tau*fold

for n in range(0,timesteps):
    #fnew = f(n*tau,p[n])
    #print "time step index", n
    #print "fnew", fnew
    #p[n+1] = p[n] + tau*(-fold/2. + 3.*fnew/2.)
    #fold = fnew
    ptemp = p[n] + tau*f(n*tau,p[n])/2.
    p[n+1] = p[n] + tau*f(tau*n+tau/2., ptemp)

print("task b) numerical scheme output:", p)

plt.plot(time,p,'go',time,p,linewidth=3)
plt.show()