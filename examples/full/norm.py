from math import sqrt

### Inputs

N0=2e19
dx=2.773e-5
B=0.55
E0=1e6
w_coeff=1.62
T=57.6 ## electron volts 

### Constants

c=299792458
## Electron charge
q_e=1.60217662e-19
## Electron mass
m_e=9.10938356e-31
## Vacuum permittivity  
eps0 = 8.85418782e-12
## Vacuum permeability
mu0 = 1.25663706e-6
##
kb=1.38064852e-23


w=w_coeff*(B*q_e/m_e)
v=sqrt(2*T*q_e / (3*m_e))

## Normalized inputs
v_std=v/c
E0_n=E0/(c**2)*sqrt(dx**3/(m_e*mu0))
B_n=B/c*sqrt(dx**3/(mu0*m_e))
w_n=w*dx/c
q_dens=q_e*N0*(dx**3)*sqrt(mu0/(m_e*dx))
m_dens=m_e*N0*(dx**3)/m_e


print("Source E-field strength",E0_n)
print("Static B-field strength",B_n)
print("Angular frequency ", w_n)
print("Standard deviation for velocity",v_std)
print("Mass density per cell",q_dens)
print("Charge density per cell",m_dens)
