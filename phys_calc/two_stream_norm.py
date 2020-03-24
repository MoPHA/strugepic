from math import sqrt
## Inputs in SI units
## Length in X direction
## Number of cells in X direction

mag=15
## Electric field strength
E=0.0
## Temp in ev
T=1e3
## Particle density
n=1e12
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




## 

## Derived values
## s_ means simulation value
# Debye
dx = sqrt(eps0*kb*T/(n*q_e*q_e))
dy = dx
dz = dx
epc = n*dx**3
w_pe=sqrt(n*q_e**2/(m_e*eps0))

## Magnetic field strength
B=mag*w_pe*m_e/q_e

v = sqrt(2*kb*T/m_e)
r = v*m_e/(q_e*B)
vd = 2.5*v

print("Physical parameters")
print("Debye length: ",dx)
print("Gyroradius: ",r)
print("Ratio: ",r/dx)
print("Plasma freq: ",w_pe)
print("Cyclotron freq:",q_e*B/m_e)
print("Ratio: ",(q_e*B/m_e)/w_pe)
print("Magnetic field strength: ",B)
print("Speed of light: ",c)
print("dx*N/(vd/w_pe)",dx*2000/(vd/w_pe))

print("\nSimulation paramters: ")
s_dx=1
s_B=B*(1/sqrt(mu0*m_e))*sqrt(dx**3)/c
s_E=E/(c**2*sqrt(mu0*m_e/dx**3))
s_v=v/c
s_vd=s_v*2.5;
s_q=q_e/(sqrt(m_e*dx/mu0))*epc
s_m=m_e/m_e*epc
s_wpe=w_pe*dx/c
print("Drift velocity: ", s_vd)
print("Thermal velocity: ", s_v)
print("Charge per cell: ",s_q)
print("Mass per cell: ", s_m)
print("Magnetic field: ", s_B)
print("Steps per period", 2/s_wpe)
print("dt*wp",0.5*s_wpe)
print("Gyroradius",s_v*s_m/(s_q*s_B))
