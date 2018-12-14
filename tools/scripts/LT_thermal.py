#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, complex_ode, ode

##############################################################################
# Reference:                                                                 #
# [1] Charles F. Gammie (1996). "Linear Theory of Magnetized, Viscous, Self- #
#     Gravitating Gas Disks"                                                 #
#     The Astrophysical Journal 462: 725-731                                 #
# [2] Charles F. Gammie (2001). "Nonlinear Outcome of Gravitational          #
#     Instability in Cooling, Gaseous Disks"                                 #
#     The Astrophysical Journal 553: 174-183                                 #
##############################################################################

Omega = 1.0                                       # angular velocity         #
TSIM = 10./Omega                                  # simulation time          #
gamma = 1.4                                       # heat capacity ratio      #
G = 6.6742e-11                                    # gravitational constant   #
Sigma0 = 1.0                                      # background density       #
L = 40.*G*Sigma0/Omega**2                         # field size               #
P0 = np.pi**2.*G**2.*Sigma0**3./(gamma*Omega**2.) # fullfills toomre's crit. #
kx0 = -2.*(2.*np.pi/L)                            # wave number              #
ky0 = 2*np.pi/L                                   # wave number              #

dsigma0 = Sigma0*5e-4                             # init. dens. pert.        #
dvx0 = 0.0                                        # init. vel. pert.         #
dvy0 = 0.0                                        # init. vel. pert.         #
dp0 = 0.0                                         # init. press. pert.       #

# wave vector
def kx(t):
  global kx0,ky0,Omega
  return kx0 + 3./2.*Omega*ky0*t
def k2(t):
  global ky0
  return ky0**(2.) + kx(t)**(2.)

# initial conditions
y0 = [dsigma0,dvx0,dvy0,dp0]
t0 = 0.0
dt = 0.001
t = []
dSIGMA = []
dVX = []
dVY = []

# definition of first order differential equation system
def ydash(t,y):
  global Sigma0, P0, Omega, gamma, G, ky0

  Sigma = y[0]                                    # density perturbation     #
  vx = y[1]                                       # xvel. perturbation       #
  vy = y[2]                                       # yvel. perturbation       #
  P = y[3]                                        # press. perturbation      #

  f0 = complex(0,-Sigma0*(kx(t)*vx + ky0*vy))
  f1 = complex(2.0*Omega*vy, kx(t)*Sigma*2.*np.pi*G/np.sqrt(k2(t))-kx(t)*P/Sigma0)
  f2 = complex(-0.5*Omega*vx,ky0  *Sigma*2.*np.pi*G/np.sqrt(k2(t))-ky0  *P/Sigma0)
  f3 = complex(0,-gamma*P0*(kx(t)*vx + ky0*vy))
  return [f0,f1,f2,f3]

# integration
r = ode(ydash).set_integrator('zvode')
r.set_initial_value(y0,t0)
while r.successful() and r.t < TSIM:
  r.integrate(r.t+dt)
  print(r.y[0].real)
  t.append(r.t)
  dSIGMA.append(r.y[0])
  dVX.append(r.y[1])
  dVY.append(r.y[2])


dSIGMA = np.array(dSIGMA)
t = np.array(t)

# plotting
plt.figure()

## command line output
#for i, k in enumerate(y[:,0]):
#  print y[i,0]

plt.plot(t*Omega,dSIGMA/Sigma0,label='Density Perturbation')
plt.xlabel('time')
plt.ylabel('amplitude')
plt.title('Linear Theory')
plt.legend(loc=0)
plt.grid(True)

plt.show()
