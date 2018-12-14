#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.signal import argrelextrema

##############################################################################
# Reference:                                                                 #
# [1] Charles F. Gammie (1996). "Linear Theory of Magnetized, Viscous, Self- #
#     Gravitating Gas Disks"                                                 #
#     The Astrophysical Journal 462: 725-731                                 #
# [2] Charles F. Gammie (2001). "Nonlinear Outcome of Gravitational          #
#     Instability in Cooling, Gaseous Disks"                                 #
#     The Astrophysical Journal 553: 174-183                                 #
# [3] Sijme-Jan Paardekooper (2012). "Numerical convergence in self-         #
#     gravitating shearing sheet simulations and the stochastic nature of    #
#     disc fragmentation"                                                    #
#     Mon. Not. R. Astron. Soc. 421, 3286-3299                               #
##############################################################################

Omega = 1.0                                       # angular velocity         #
TSIM = 10. / Omega                                # simulation time          #
G = 1.0	                                          # gravitational constant   #
# background
Sigma0 = 1.0 / 40.
dSigma = Sigma0 * 5e-4
L = 1.                                            # field size               #
kx0 = -4. * np.pi                                 # wave number              #
ky0 = 2. * np.pi                                  # wave number              #
Ac = -0.75 * Omega                                # Oort's constant          #
# Oort's constant          #
Bc = (Omega + Ac)
kappa = Omega                                     # epicyclic frequency      #
# speed of sound           #
cs0 = np.pi * G * Sigma0 / kappa
# ATTENTION: Depends on Sigma0!!
dxi0 = -2.5e-4/Sigma0                             # initial xi-perturbation  #

# background vorticity     #
xi0 = 2. * Bc / Sigma0
dtSigma0 = 0.0                                    # init. vel. dens. pert.   #
nu = 0.0


# wavenumbers
def kx(t):
    global kx0, ky0, Ac
    return kx0 - 2 * Ac * ky0 * t


def k2(t):
    global ky0
    return ky0**(2.) + kx(t)**(2.)

# functions in Hunter's equation


def a(t):
    global Ac, ky0, nu
    return 4. * Ac * kx(t) * ky0 / k2(t) + k2(t) * 4. * nu / 3.


def b(t):
    global kappa, G, Sigma0, cs0, Ac, Bc, ky0
    return (kappa**(2.) - 2. * np.pi * G * Sigma0 * np.sqrt(k2(t)) +
            8. * Ac * Bc * ky0**2. / k2(t))


def c(t):
    global Sigma0, Omega, Ac, ky0
    return Sigma0**(2.) * (2. * Omega + 4. * Ac * ky0**(2.) / (k2(t)))


def c2(t):
    global Sigma0, cs0
    return 2. * cs0 * k2(t)

# functions in vorticity equation


def d(t):
    global nu
    return nu * k2(t)


def e(t):
    global xi0, nu, Sigma0
    return nu * k2(t) * xi0 / Sigma0

# function in energy/soundspeed equation


def f(t):
    global gamma, cs0
    return (gamma - 1.) / 2.

# definition of first order differential equation system


def ydash(y, t):
    global Sigma0, cs0, gamma, P0, dxi0
    # density perturbation     #
    sig = y[0]
    # derivation of dens. pert.#
    dsig = y[1]

    f0 = dsig
#  f1 = - dsig*a(t) - sig*(b(t) + gamma*P0*k2(t)) - xi*c(t) + c1*gamma*P0
    f1 = - dsig * a(t) - sig * (b(t) + cs0**(2.) * k2(t)) - dxi0 * c(t)
#  f2 = - xi*d(t) - sig*e(t)
#  return [f0,f1,f2]
    return [f0, f1]


plt.figure()

# initial conditions
# y0 = [5e-4*Sigma0,dtSigma0,dxi0]
y0 = [5e-4 * Sigma0, dtSigma0]
t = np.linspace(0, TSIM, 10000)
# solve ode
y = odeint(ydash, y0, t)

# command line output
np.savetxt('lt_isothermal.dat', y[:, :]/Sigma0)
print("Results written to lt_isothermal.dat")
plt.plot(t * Omega, y[:, 0] / Sigma0, label='Density')

# Give out local maxima and minima (time and value)
print("Maxima:", "time:", t[argrelextrema(y[:, 0] / Sigma0, np.greater)] *
      Omega, y[argrelextrema(y[:, 0] / Sigma0, np.greater),0]/Sigma0)
print("Minima:", "time:", t[argrelextrema(y[:, 0] / Sigma0, np.less)]*Omega,
      y[argrelextrema(y[:, 0] / Sigma0, np.less), 0]/Sigma0)

plt.xlabel('time')
plt.ylabel('amplitude')
plt.title('Linear Theory')
plt.legend(loc=0)
plt.grid(True)

plt.show()
