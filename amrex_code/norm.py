## Inputs in SI units
## Length in X direction
L=1.5e-4
## Number of cells in X direction
N=64
## Magnetic field strength
B=0.55
## Particle velocity
v=6e6
## Particle density
n=1
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

## Derived values
## s_ means simulation value

dx = L/N
dy = dx
dz = dx
epc = n*dx**3
qc_e= q_e*epc
mc_e= m_e*epc


s_dx=1
s_L=N
s_dt=dt*c/dx
s_B=B/(c*dx*mu0)
s_v=v/c
s_q=qc_e*c/(dx**3)
s_m=mc_e*c**3*eps0/(dx**5)

print("Real radius: " + str(v*m_e/(q_e*B)))
print("Real omega: " +str(q_e*B/m_e))
print("Time around: " + str( (2*3.1415962)/(q_e*B/m_e) ))
print("Time steps around: " + str(int( (2*3.1415962)/(q_e*B/m_e*dt) )))
print("Scaled radius: " +str(s_v*s_m/(s_q*s_B)))


print("Input parameters: ")
print("Simulation length X: ",s_L)
print("Simulation length Y: ",s_L)
print("Magnetic field: " ,s_B)
print("Initial velocity: " ,s_v)
print("Charge: " ,s_q)
print("Mass: " ,s_m)

