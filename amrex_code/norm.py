from math import sqrt
## Inputs in SI units
## Length in X direction
L=4e-2
## Number of cells in X direction
N=1800
## Magnetic field strength
B=0.55
## Electric field strength
E=1e6
## Particle velocity
v=6e6
## Temp in ev
T= 57.6
## Particle density
n=2e19
## Time step
dt=9.24e-14
## Number of time steps
Nt=1000
## speed of light
c=299792458
## Electron charge
q_e=1.60217662e-19
## Electron mass
m_e=9.10938356e-31
## Vacuum permittivity  
eps0 = 8.85418782e-12
##
mu0 = 1.25663706e-6

## 
E_omega=1.62*9.67e10

## Derived values
## s_ means simulation value

dx = L/N
dy = dx
dz = dx
epc = n*dx**3

s_dt=dt*c/dx
s_dx=1
s_L=N
s_B=B/(c*dx*mu0)
s_E=E*eps0/dx
s_v=sqrt(T/(m_e*c**2*6.242e+18))

s_l = (c/E_omega)/dx
s_omega = dx*E_omega/c
s_q=q_e/(dx**3)*epc
s_m=m_e/(mu0*dx**5)*epc


print("\nInput parameters: ")
print("Simulation length X: ",s_L)
print("Magnetic field: " ,s_B)
print("Electric field: " ,s_E)
print("Wave freq: " ,s_omega)
print("Cyclotron freq:", (s_B*(q_e/(dx**3) )/(m_e/(mu0*dx**5))  ))
print("Freq ratio: ", (dx*E_omega/c)/(s_B*(q_e/(dx**3) )/(m_e/(mu0*dx**5))  ))
print("Initial velocity (Maxwell stand.dev): " ,s_v)
print("Charge per cell: " ,s_q)
print("Mass per cell: " ,s_m)



