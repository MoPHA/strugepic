import matplotlib.pyplot as plt
import numpy as np




arr=np.load('data.npy')
plt.imshow(arr[:,:],interpolation='gaussian', cmap='BrBG',origin='lower',aspect='auto')
plt.show()
