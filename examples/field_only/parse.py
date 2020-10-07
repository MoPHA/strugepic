import yt
from os import path
import numpy as np
import matplotlib.pyplot as plt
from numpy import save

steps=range(0,4000,10)
evolution=np.zeros((len(steps),1000))


i=0
for step in steps:
    E_filename = "Field_Only_Data_out/plt_E" +str(step)
    dsE=yt.load(E_filename)
    E=dsE.r[:]
    Ey=np.array(E['E_y'])
    x =np.array( E['x'])
    y =np.array( E['y'])
    z =np.array( E['z'])
    data=np.zeros((1000*10*1,4))
    data[:,0]=y-0.5
    data[:,1]=x-0.5
    data[:,2]=z-0.5
    data[:,3]=Ey
    data = data[data[:,2].argsort()] # First sort doesn't need to be stable.
    data = data[data[:,1].argsort(kind='mergesort')] 
    data = data[data[:,0].argsort(kind='mergesort')]
    E_y=data[:,3].reshape(10,1000)
    E_y_l=E_y[5,:]
    evolution[i,:]=E_y_l
    print(i) 
    i=i+1
save('data.npy', evolution)
