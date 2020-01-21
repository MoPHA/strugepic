import yt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

X = np.linspace(-1,1, 256)
Y = np.linspace(-1,1, 256)
X, Y = np.meshgrid(X, Y)

#for step in range(0,70):
#    filename = "plt_E"+str(step)
#    ds=yt.load(filename)
#    mag=ds.r[:,:,0] 
#    arr =mag["E_y"].reshape(256,256)
#    #fig = plt.figure()
#    #ax = fig.gca(projection='3d')
#
#    plt.imshow(arr, cmap='summer', interpolation='nearest')
#    #ax.plot_surface(X,Y,arr,cmap=cm.coolwarm)
#    #ax.set_ylim([-1,1])
#    #ax.set_xlim([-1,1])
#
#    print(step)
#    plt.savefig("E_y_field"+str(step).zfill(3)) 
#    plt.close()

for step in range(0,1000,1):
    B_filename = "plt_B" + str(step)
    E_filename = "plt_E" +str(step)
    dsB=yt.load(B_filename)
    dsE=yt.load(E_filename)
    B=dsB.r[[0,2,2]:[512,2,2]]
    E=dsE.r[[0,2,2]:[512,2,2]]
    #B=dsB.r[:,:,:]
    #E=dsE.r[:,:,:]
    Bx = B["B_x"]
    By = B["B_y"]
    Bz = B["B_z"]
    Ex = E["E_x"]
    Ey = E["E_y"]
    Ez = E["E_z"]


    fig, ax = plt.subplots()
    ax.plot(Bz)
    ax.plot(Ey)
    ax.set_ylim([-1,1])
 #   ax.set_title(str(val[i-1]))
    plt.savefig("B_pic"+str(step).zfill(3)) 
    plt.close()
 #   fig, ax = plt.subplots()
 #   ax.plot(Ex)
 #   ax.plot(Ey)
 #   ax.plot(Ez)
 #   ax.set_ylim([-1,1])
 #   ax.set_title(str(val[i-1]))
 #   plt.savefig("E_pic"+str(step).zfill(3)) 
 #   plt.close()


