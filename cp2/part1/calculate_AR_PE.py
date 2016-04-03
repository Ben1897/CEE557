import csv
import numpy as np
import matplotlib.pyplot as plt

# parameters
dx = 0.05    # space resolution [m]
dt = 9     # time resolution [day]
D  = 2.5e-5  # dispersion coefficient [m2/day]
u  = 0.005   # advection speed [m/day]

SD = D*dt/(dx**2)  # script D
SU = u*dt/dx       # script U

wldx = np.arange(2., 100., 2)  # the ratio of the wavelength over the space resolution

# function for writing the AR and PE results
def writeARPE(filename, wldx, AR, PE):
    with open(filename, 'w') as csvfile:
        fieldnames = ['wavelength / dx', 'AR', 'PE']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for i in xrange(len(wldx)):
            writer.writerow({'wavelength / dx':      "%.1f" % wldx[i],
                             'AR':                   "%.7f" % AR[i],
                             "PE":                   "%.7F" % PE[i]})

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
def AF_fdm_CN_central(wldx):
    values = []
    for ele in wldx:
        value = complex(1-SD*(1.-np.cos(2.*np.pi/ele)), -SU/2*np.sin(2*np.pi/ele)) / \
                complex(1+SD*(1.-np.cos(2.*np.pi/ele)), SU/2*np.sin(2*np.pi/ele))
        values.append(value)
    return np.array(values)
AF_fdm_CN_central_mod = lambda wldx: np.absolute( AF_fdm_CN_central(wldx) )
PS_fdm_CN_central = lambda wldx: np.arcsin(np.imag(AF_fdm_CN_central(wldx)) /
                                           AF_fdm_CN_central_mod(wldx))
AR_fdm_CN_central = lambda wldx: (AF_fdm_CN_central_mod(wldx) /
                                  AF_anal_mod(wldx)) ** Nn(wldx)
PE_fdm_CN_central = lambda wldx: -2*np.pi - Nn(wldx) * PS_fdm_CN_central(wldx)
# upstream
def AF_fdm_CN_up(wldx):
    values = []
    for ele in wldx:
        value = complex(1-(SD+SU/2.)*(1.-np.cos(2.*np.pi/ele)), -SU/2*np.sin(2*np.pi/ele)) / \
                complex(1+(SD+SU/2.)*(1.-np.cos(2.*np.pi/ele)), SU/2*np.sin(2*np.pi/ele))
        values.append(value)
    return np.array(values)
AF_fdm_CN_up_mod = lambda wldx: np.absolute( AF_fdm_CN_up(wldx) )
PS_fdm_CN_up = lambda wldx: np.arcsin(np.imag(AF_fdm_CN_up(wldx)) /
                                      AF_fdm_CN_up_mod(wldx))
AR_fdm_CN_up = lambda wldx: (AF_fdm_CN_up_mod(wldx) /
                             AF_anal_mod(wldx)) ** Nn(wldx)
PE_fdm_CN_up = lambda wldx: -2*np.pi - Nn(wldx) * PS_fdm_CN_up(wldx)

# FEM

# plot & write the results to the table
# AR & PE plots of upstream and central differentce methods by FDM-explicit (Part A - 1(a))
AR_fdm_explicit_central_values = AR_fdm_explicit_central(wldx)
PE_fdm_explicit_central_values = PE_fdm_explicit_central(wldx)
AR_fdm_explicit_up_values = AR_fdm_explicit_up(wldx)
PE_fdm_explicit_up_values = PE_fdm_explicit_up(wldx)
fig, (ax11, ax12) = plt.subplots(2, 1)
ax11.plot(wldx, AR_fdm_explicit_central_values, '.--', label="central")
ax11.plot(wldx, AR_fdm_explicit_up_values, '-', label="upstream")
ax11.set_title('FDM - explicit')
ax11.set_ylabel('A.R.')
ax11.legend(loc='upper right')
ax12.plot(wldx, PE_fdm_explicit_central_values, '.--', label="central")
ax12.plot(wldx, PE_fdm_explicit_up_values, '-', label="upstream")
ax12.set_xlabel('wavelength/dx')
ax12.set_ylabel('P.E.')
ax12.legend(loc='lower right')
plt.savefig('FDM_explicit_AR_PE.png')
writeARPE('FDM_explicit_central_AR_PE.csv', wldx, AR_fdm_explicit_central_values, PE_fdm_explicit_central_values)
writeARPE('FDM_explicit_up_AR_PE.csv', wldx, AR_fdm_explicit_up_values, PE_fdm_explicit_up_values)

# AR & PE plots of upstream and central differentce methods by FDM-Crank Nicolson (Part A - 2(a))
AR_fdm_CN_central_values = AR_fdm_CN_central(wldx)
PE_fdm_CN_central_values = PE_fdm_CN_central(wldx)
AR_fdm_CN_up_values = AR_fdm_CN_up(wldx)
PE_fdm_CN_up_values = PE_fdm_CN_up(wldx)
fig, (ax21, ax22) = plt.subplots(2, 1)
ax21.plot(wldx, AR_fdm_CN_central_values, '.--', label="central")
ax21.plot(wldx, AR_fdm_CN_up_values, '-', label="upstream")
ax21.set_title('FDM - Crank Nicolson')
ax21.set_ylabel('A.R.')
ax21.legend(loc='upper right')
ax22.plot(wldx, PE_fdm_CN_central_values, '.--', label="central")
ax22.plot(wldx, PE_fdm_CN_up_values, '-', label="upstream")
ax22.set_xlabel('wavelength/dx')
ax22.set_ylabel('P.E.')
ax22.legend(loc='lower right')
plt.savefig('FDM_CN_AR_PE.png')
writeARPE('FDM_CN_central_AR_PE.csv', wldx, AR_fdm_CN_central_values, PE_fdm_CN_central_values)
writeARPE('FDM_CN_up_AR_PE.csv', wldx, AR_fdm_CN_up_values, PE_fdm_CN_up_values)

fig, (ax51, ax52) = plt.subplots(2, 1)
ax51.plot(wldx, AF_fdm_explicit_central_mod(wldx), '.--', label="central")
ax51.plot(wldx, AF_fdm_explicit_up_mod(wldx), '-', label="upstream")
ax51.plot(wldx, AF_fdm_CN_central_mod(wldx), '.-.', label="CN-central")
ax51.plot(wldx, AF_fdm_CN_up_mod(wldx), '*-', label="CN-upstream")
ax51.plot(wldx, AF_anal_mod(wldx), 'o-', label="analytical")
ax51.legend(loc='lower right')
ax52.plot(wldx, PS_fdm_explicit_central(wldx), '.--', label="central")
ax52.plot(wldx, PS_fdm_explicit_up(wldx), '-', label="upstream")
ax52.plot(wldx, PS_fdm_CN_central(wldx), '.-.', label="CN-central")
ax52.plot(wldx, PS_fdm_CN_up(wldx), '*-', label="CN-upstream")
ax52.plot(wldx, PS_anal(wldx), 'o-', label="analytical")
ax52.legend(loc='lower right')

plt.show()
# plt.close('all')
