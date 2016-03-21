import numpy as np
import matplotlib.pyplot as plt

# parameters
dx = 0.05    # space resolution [m]
dt = 1.0     # time resolution [day]
D  = 2.5e-5  # dispersion coefficient [m2/day]
u  = 0.005   # advection speed [m/day]

SD = D*dt/(dx**2)  # script D
SU = u*dt/dx       # script U

wldx = np.arange(2, 60, 2)  # the ratio of the wavelength over the space resolution

# functions for calculating
# (1) the modular of the amplification factor
# (2) the phase shift

Nn = lambda wldx: wldx/SU

# analytical solution
AF_anal_mod = lambda wldx: np.exp(-4*np.pi**2*SD*(1/wldx)**2)  # amplification factor
PS_anal = lambda wldx: -2*np.pi*SU*(1/wldx)  # phase shift

# FDM - explicit
# central
AF_fdm_explicit_central_mod = lambda wldx: np.sqrt((1-2*SD*(1-np.cos(2*np.pi/wldx)))**2 +
                                                   (SU*np.sin(2*np.pi/wldx))**2)
PS_fdm_explicit_central = lambda wldx: np.arcsin(-SU*np.sin(2*np.pi/wldx) /
                                                 AF_fdm_explicit_central_mod(wldx))
AR_fdm_explicit_central = lambda wldx: (AF_fdm_explicit_central_mod(wldx) /
                                        AF_anal_mod(wldx)) ** Nn(wldx)
PE_fdm_explicit_central = lambda wldx: -2*np.pi - Nn(wldx) * PS_fdm_explicit_central(wldx)
# upstream
AF_fdm_explicit_up_mod = lambda wldx: np.sqrt((1-(2*SD+SU)*(1-np.cos(2*np.pi/wldx)))**2 +
                                              (SU*np.sin(2*np.pi/wldx))**2)
PS_fdm_explicit_up = lambda wldx: np.arcsin(-SU*np.sin(2*np.pi/wldx) /
                                            AF_fdm_explicit_up_mod(wldx))
AR_fdm_explicit_up = lambda wldx: (AF_fdm_explicit_up_mod(wldx) /
                                   AF_anal_mod(wldx)) ** Nn(wldx)
PE_fdm_explicit_up = lambda wldx: -2*np.pi - Nn(wldx) * PS_fdm_explicit_up(wldx)

# FDM - Crank-Nilcoson
# central
# upstream

# FEM

# plot
# AR & PE plots of upstream and central differentce methods by FDM (Part A - (a))
AR_fdm_explicit_central_values = AR_fdm_explicit_central(wldx)
PE_fdm_explicit_central_values = PE_fdm_explicit_central(wldx)
AR_fdm_explicit_up_values = AR_fdm_explicit_up(wldx)
PE_fdm_explicit_up_values = PE_fdm_explicit_up(wldx)
fig, ax = plt.subplots()
ax.plot(wldx, AR_fdm_explicit_central_values, '.--', label="central")
ax.plot(wldx, AR_fdm_explicit_up_values, '-', label="upstream")
ax.legend(loc='lower right')

fig, ax2 = plt.subplots()
ax2.plot(wldx, PE_fdm_explicit_central_values, '.--', label="central")
ax2.plot(wldx, PE_fdm_explicit_up_values, '-', label="upstream")
ax2.legend(loc='lower right')

plt.show()
