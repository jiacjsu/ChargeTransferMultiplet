import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

pols = ['xinyout', 'xinzout', 'yinxout', 'yinzout', 'zinxout', 'zinyout']

fig = plt.figure(figsize=(12,6))
mydir = '/Users/cjjia/Documents/Research/MatthiasRIXS/LaNiO2_RIXS/'
data_nsp = np.loadtxt(mydir + 'LaNiO2.RIXS-nsf.dat')
data_sp = np.loadtxt(mydir + 'LaNiO2.RIXS-sf.dat')

divX = 21 #incoming energy
divY = 101 # energy loss
Z1 = data_nsp[:,2].reshape(divX,divY) 
Z2 = data_sp[:,2].reshape(divX,divY)

im = plt.imshow(np.transpose(Z1 + Z2), interpolation='bilinear', extent=[851, 853,0, 5],
                origin='lower',vmax=(Z1 + Z2).max()/1.0, vmin=0, cmap = 'coolwarm')

plt.colorbar(orientation ='vertical')
plt.show()





