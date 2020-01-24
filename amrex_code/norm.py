## Inputs in SI units
## Length in X direction
L=4e-2
## Number of cells in X direction
N=1800
## Magnetic field strength
B=0.5
## Electric field strength
E=1e6
## Particle velocity
v=6e6
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
s_v=v/c

s_l = (c/E_omega)/dx
s_omega = dx*E_omega/c
s_q=q_e*c/(dx**3)
s_m=m_e*c**3*eps0/(dx**5)


print("\nInput parameters: ")
print("Simulation length X: ",s_L)
print("Simulation length Y: ",s_L)
print("Particle per cell: " , epc)
print("Magnetic field: " ,s_B)
print("Electric field: " ,s_E)
print("Wave freq: " ,s_omega)
print("Initial velocity: " ,s_v)
print("Charge: " ,s_q)
print("Mass: " ,s_m)

