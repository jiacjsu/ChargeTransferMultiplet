# concerns: right shift matters for signed or unsigned int? 

# input parameters for GenHspCoreFull:
#   N_d, N_ligand, N_core, nup, ndn: all are integers
#   Hsp = np.zeros((10000000), dtype = np.uint32) 
#   GenHsp2p1Hole(8,11,11,Hsp) for Ferricyanide


# GenHsp2p1Hole speed still too slow for Kristjan's problem Ferricyanide

import numpy as np
from scipy.special import comb
import ctypes

def sumEveryBit(tempa):
    tempsum=0
    tempa1=tempa
    while tempa1:
        tempsum += tempa1 % 2
        tempa1 = tempa1 >> 1
    return tempsum


def GenHspCoreFull(N_ligand, nup, ndn, Hsp):

    N_core = 3
    N_d = 5    
    N = N_d + N_ligand + N_core
    ntemp = 10000
    istatup = np.zeros((ntemp), dtype = np.uint64)
    istatdn = np.zeros((ntemp), dtype = np.uint64)
    
    if N_core == 3:
        flagCoreFull = 7 * (1 << (N-N_core)) 
    
    upsize = 0
    for istat in range(1 << N):
        if(sumEveryBit(istat) == nup and (istat & flagCoreFull) == flagCoreFull):
             istatup[upsize] = istat
             upsize += 1
    
    dnsize = 0
    for istat in range(1 << N):
        if(sumEveryBit(istat) == ndn and (istat & flagCoreFull) == flagCoreFull): 
            istatdn[dnsize] = istat
            dnsize += 1
    print '        GenHspCoreFull    upsize, dnsize: ', upsize, dnsize

    Hsize = 0
    for ii in range(upsize):
        for jj in range(dnsize):
            #print 'statup, statdn: ', istatup[ii], istatdn[jj]
            Hsp[Hsize] = istatdn[jj] + istatup[ii] * 2**N
            Hsize += 1
    #Hsp = Hsp[:Hsize]
    print ''
    return Hsize



#!************************************************

# dipole excitation is always for spin-up electrons.
# dipole de-excitation could be spin-up or spin-down electrons.

def GenHsp2p1Hole(N_ligand, nup, ndn, Hsp):

    #N_core = 3
    #N_d = 5
    #N = N_d + N_ligand + N_core
    #Hsize_i = 6 * comb(2*N-6, nup+ndn-5, exact=True)

    lib = ctypes.cdll['./GenHspCore.so']
    GenHsp2p1HoleCore = lib['GenHsp2p1HoleCore']

    Hsize = GenHsp2p1HoleCore(N_ligand, nup, ndn, ctypes.c_void_p(Hsp.ctypes.data))

    return Hsize

    
   # #print 'N ', N
   # ntemp = 10000
   # istatup = np.zeros((ntemp), dtype = np.uint32)
   # istatdn_ncore2 = np.zeros((ntemp), dtype = np.uint32)
   # istatdn_ncore3 = np.zeros((ntemp), dtype = np.uint32)
   # 
   # if N_core == 3:
   #     flagCoreFull = 7 * (1 << (N-N_core)) 
   # #print 'flagCoreFull ', bin(flagCoreFull)
   #     
   # upsize = 0
   # for sup in range(1 << N): 
   #     #if sup % 1024 == 0: print sup
   #     if (sumEveryBit(sup & flagCoreFull) == 2 and (sumEveryBit(sup) == nup)) or \
   #        (sumEveryBit(sup & flagCoreFull) == 3 and (sumEveryBit(sup) == nup + 1)):
   #         istatup[upsize] = sup
   #         upsize += 1
   # print '        GenHsp2p1Hole    upsize: ', upsize
   # istatup = istatup[:upsize] 
   # #print 'istatup[0] ', bin(istatup[0])
   # 
   # dnsize_ncore2 = 0
   # dnsize_ncore3 = 0
   # for sdn in range(1 << N): 
   #     #if sdn % 1024 == 0: print sdn
   #     if (sumEveryBit(sdn & flagCoreFull) == 3 and (sumEveryBit(sdn) == ndn)):
   #         istatdn_ncore3[dnsize_ncore3] = sdn
   #         dnsize_ncore3 += 1            
   #     if (sumEveryBit(sdn & flagCoreFull) == 2 and (sumEveryBit(sdn) == ndn - 1)):
   #         istatdn_ncore2[dnsize_ncore2] = sdn
   #         dnsize_ncore2 += 1
   # print '        GenHsp2p1Hole    dnsize_ncore2, dnsize_ncore2: ', \
   #     dnsize_ncore2, dnsize_ncore3
   # istatdn_ncore2 = istatdn_ncore2[:dnsize_ncore2]
   # istatdn_ncore3 = istatdn_ncore3[:dnsize_ncore3]
   # #print 'istatdn_ncore2[0] ', bin(istatdn_ncore2[0])
   # #print 'istatdn_ncore3[0] ', bin(istatdn_ncore3[0])
   # 
   # Hsize = 0
   # for ii in range(upsize):
   #     #if ii % 64 == 0: print ii, Hsize
   #     sup = istatup[ii]
   #     #print 'sumEveryBit ', sumEveryBit(sup)
   #     if sumEveryBit(sup & flagCoreFull) == 3:
   #         Hsp[Hsize: Hsize+dnsize_ncore2] += sup * 2**N
   #         Hsp[Hsize: Hsize+dnsize_ncore2] += istatdn_ncore2[:dnsize_ncore2] 
   #         Hsize += dnsize_ncore2
   #     if sumEveryBit(sup & flagCoreFull) == 2:
   #         Hsp[Hsize: Hsize+dnsize_ncore3] += sup * 2**N
   #         Hsp[Hsize: Hsize+dnsize_ncore3] += istatdn_ncore3[:dnsize_ncore3] 
   #         Hsize += dnsize_ncore3
   #         '''
   #     for jj in range(dnsize):
   #         sup = istatup[ii]
   #         sdn = istatdn[jj]
   #         #print 'sup, sdn, type(sup), sumEveryBit(sup): ', sup, sdn, type(sup), sumEveryBit(sup)
   #         if (sumEveryBit(sup)+sumEveryBit(sdn)) == (nup+ndn):   # and \
   #         #   (sumEveryBit(sup & flagCoreFull)+sumEveryBit(sdn & flagCoreFull)) == 5:
   #             Hsp[Hsize] = sdn + sup * 2**N
   #             Hsize += 1
   #         '''
   # #Hsp = Hsp[:Hsize]
   # #print 'Hsp[0]', bin(Hsp[0])
   # print ''
   # return Hsize
