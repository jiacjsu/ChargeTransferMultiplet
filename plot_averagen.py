import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

pols = ['xinyout', 'xinzout', 'yinxout', 'yinzout', 'zinxout', 'zinyout']

fig = plt.figure(figsize=(16,8))



data_n = np.loadtxt('/Users/cjjia/Documents/Research/MatthiasRIXS/NiO_ddonly_checkfinalstates/NiO.averagen-nsf.dat')
width = 0.01

divX = 51 #incoming energy
divY = 101 # energy loss

plt.subplot(4,1,1)
plt.bar(data_n[:,0], data_n[:,1], width, label = 'dz2 up', color = 'r')
plt.legend()

plt.subplot(4,1,2)
plt.bar(data_n[:,0], data_n[:,2], width, label = 'dx2-y2 up', color = 'g')
plt.legend()

plt.subplot(4,1,3)
plt.bar(data_n[:,0], data_n[:, 1+5+3+3], width, label = 'dz2 dn', color = 'b')
plt.legend()

plt.subplot(4,1,4)
plt.bar(data_n[:,0], data_n[:, 2+5+3+3], width, label = 'dx2-y2 dn', color = 'm')
plt.legend()

plt.title('4up4dn + sf')
plt.legend()



plt.title('Total')

plt.show()





