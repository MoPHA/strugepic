import yt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


for step in range(0,700,5):
    filename = "plt_E" + str(step)
    ds=yt.load(filename)
    mag=ds.r[[-1,0,0]:[1,0,0]]
    fig, ax = plt.subplots()
    ax.plot(mag["E_z"])
    ax.set_ylim([-1,1])
    plt.savefig("E_z_pic"+str(step).zfill(3)) 
    plt.close()
