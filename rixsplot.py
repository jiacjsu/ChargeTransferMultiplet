import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

pols = ['xinyout', 'xinzout', 'yinxout', 'yinzout', 'zinxout', 'zinyout']

fig = plt.figure(figsize=(16,8))



data_nsp = np.loadtxt('/Users/cjjia/Documents/Research/MatthiasRIXS/src_callC_NiO_64bit_Ecorr/NiO.RIXS-nsf.dat')
data_sp = np.loadtxt('/Users/cjjia/Documents/Research/MatthiasRIXS/src_callC_NiO_64bit_Ecorr/NiO.RIXS-sf.dat')

divX = 51 #incoming energy
divY = 101 # energy loss
Z1 = data_nsp[:,2].reshape(divX,divY) 
Z2 = data_sp[:,2].reshape(divX,divY)

ax1 = plt.subplot(2,4,1)
im = plt.imshow(np.transpose(Z1), interpolation='bilinear', extent=[855, 860, 0, 5],
                origin='lower',vmax=Z1.max()/2.0, vmin=0, cmap = 'coolwarm')
plt.colorbar(orientation ='vertical')
plt.title('spin down + nsf')

ax2 = plt.subplot(2,4,2)
im = plt.imshow(np.transpose(Z2), interpolation='bilinear', extent=[855, 860, 0, 5],
                origin='lower',vmax=Z2.max()/2.0, vmin=0, cmap = 'coolwarm')
plt.colorbar(orientation ='vertical')
plt.title('spin down + sf')


# total
ax5 = plt.subplot(1,2,2)
Ztol = Z1 + Z2 
im = plt.imshow(np.transpose(Z1), interpolation='bilinear', extent=[855, 860, 0, 5],
                origin='lower',vmax=Ztol.max()/2.0, vmin=0, cmap = 'coolwarm')
plt.colorbar(orientation ='vertical')
plt.title('Total')

plt.show()





