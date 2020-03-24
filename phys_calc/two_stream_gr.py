import numpy as np
import matplotlib.pyplot as plt
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
kb=1.38064852e-23



## Temperature
T=1e3
## Particle density
n=1e12

w_pe=np.sqrt(n*q_e**2/(m_e*eps0))
v_drift=np.sqrt(2*kb*T/m_e)*2.5
print(v_drift)
k=(np.linspace(-2,2,1200)+0j)*w_pe/v_drift
x= k*v_drift/w_pe

w1 = w_pe*np.sqrt((1+2*x**2+np.sqrt(1+8*x**2))/2)/w_pe
w2 = w_pe*np.sqrt((1+2*x**2-np.sqrt(1+8*x**2))/2)/w_pe
w3 = -w_pe*np.sqrt((1+2*x**2-np.sqrt(1+8*x**2))/2)/w_pe
w4 = -w_pe*np.sqrt((1+2*x**2+np.sqrt(1+8*x**2))/2)/w_pe

#w1=  np.sqrt(k**2*v_drift**2+w_pe**2+w_pe*np.sqrt(4*k**2*v_drift**2+w_pe**2))/w_pe 
#w2=  np.sqrt(k**2*v_drift**2+w_pe**2-w_pe*np.sqrt(4*k**2*v_drift**2+w_pe**2))/w_pe 
#w3= -np.sqrt(k**2*v_drift**2+w_pe**2-w_pe*np.sqrt(4*k**2*v_drift**2+w_pe**2))/w_pe 
#w4= -np.sqrt(k**2*v_drift**2+w_pe**2+w_pe*np.sqrt(4*k**2*v_drift**2+w_pe**2))/w_pe 
#
w2I = w2.imag
w2R = w2.real
w3I = w3.imag
w3R = w3.real
#w2R[ w2.real==0 ] = np.nan
#w3I[ w3.imag==0 ] = np.nan
#w3R[ w3.real==0 ] = np.nan
#w2I[ w2.imag==0 ] = np.nan

k_scaled=k*v_drift/w_pe
fig, ax = plt.subplots()
ax.plot(k_scaled,w1,'b--')
ax.plot(k_scaled,w2I,'r-.')
x1, y1 = [-2, 2], [0, 0]
x2, y2 = [0, 0], [-2, 2]
plt.plot(x1, y1,'k')
plt.plot(x2,y2,'k')
ax.plot(k_scaled,w2R,'b--')
ax.plot(k_scaled,w3I,'r-.')
ax.plot(k_scaled,w3R,'b--')
ax.plot(k_scaled,w4,'b--')
ax.set_aspect('equal', 'box')
ax.set(ylim=(-2,2))
left,right = ax.get_xlim()
low,high = ax.get_ylim()
plt.text(-0.2,1.8,r'$\omega$',fontsize=20)
plt.text(1.8,-0.2,r'$k$',fontsize=20)
plt.ylabel(r'$\frac{\omega}{\omega_p}$',fontsize=25,rotation=0)
plt.xlabel(r'$kv_d/\omega_p$',fontsize=20)
plt.title("Dispersion diagram for two equal opposing streams",fontsize=25)
ax.legend([r'$\omega$ real part',r'$\omega$ imaginary part' ],fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=15)
plt.show()

