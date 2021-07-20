
'''
!From the LSB to HSB
!Orbital 1=d_{3z^2-r^2}; 2=d_{x^2-y^2}
!Orbital 3: dxy; 4: dxz; 5: dyz.
!Orbital 6: px; 7: py; 8: pz;
!Orbital 9-11:  O2p x+, px,py,pz
!Orbital 12-14: O2p x-, px,py,pz
!Orbital 15-17: O2p y+, px,py,pz
!Orbital 18-20: O2p y-, px,py,pz
!Orbital 21-23: O2p z+, px,py,pz
!Orbital 24-26: O2p z-, px,py,pz

!use to integers to represent spin up and spin down
!E_0(n+1 or 5) {Fe only} = -4216.93
!E_0(n or 4) {Fe only} = -4212.65
!Delta = 2.7
!E_p = E_0(n+1) - E_0(n) - Delta = -6.98
'''

    
'''
open(unit=10101,file='input',Status='old');
read(10101,'(2I8)') nup, ndn
read(10101,'(2I8)') pp, spin
read(10101,'(2F8.2, I8)') minE,maxE,Ediv
read(10101,'(2F8.2, I8)') minE_loss,maxE_loss,Ediv_loss
read(10101,'(I8)') niter_CFE
read(10101,'(2F8.2)') epsilone_CFE, epsilone_CG
close(10101)


write(*,*) ''
write(*,*) '     nup     ndn'
write(*,'(2I8)') nup, ndn
write(*,*) ''
write(*,*) '   pp     spin'
write(*,'(2I8)') pp, spin
write(*,*) ''
write(*,*) '    maxE    minE'
write(*,'(2F8.2)') maxE,minE
write(*,*) ''
write(*,*) '    Ediv'
write(*,'(I8)') Ediv
write(*,*) ''
write(*,*) '    maxEloss    minEloss'
write(*,'(2F8.2)') maxE_loss,minE_loss
write(*,*) ''
write(*,*) '    Edivloss'
write(*,'(I8)') Ediv_loss
write(*,*) ''
write(*,*) 'niter_CFE'
write(*,'(I8)') niter_CFE
write(*,*) ''
write(*,*) 'epsilone_CFE'
write(*,'(2F8.2)') epsilone_CFE, epsilone_CG
write(*,*) ''
write(*,*) '******************************'
write(*,*) ''

!============== The Step 0: preparation ==============

startY = minE; endY = maxE; divY = Ediv-1
startX = minE_loss; endX = maxE_loss;  divX = Ediv_loss-1

write(*,*) 'startX, endX ', startX, endX
write(*,*) 'startY, endY ', startY, endY
'''

import os
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg
#import scipy as sp
import numpy as np
from scipy.special import comb

from ClassDefine import ModelParas, SpectraParas
from SetModelParas import SetXASPolarParas, SetModelParas
from GenHsp import GenHspCoreFull, GenHsp2p1Hole
from GenMatrix import GenOStatesSpinUp, GenODaggerStatesSpinUpDown, GenMatrix, averageElectron
from GenMatrix import GenODaggerStatesSpinUp, GenODaggerStatesSpinDown
from CmplxContFracExpan import CmplxContFracExpan
from BiCGS import BiCGS
import matplotlib.pyplot as plt
#from ReadFromInterface import ReadFromInterface
#from ClassDefine import ModelParas, SpectraParas

#para, specPara = ReadFromInterface()

def test_Tanabe_Sugano(para):

    tenDq = np.arange(0,3.01,0.01)
    #tenDq = np.arange(0,0.01,0.01)
    data = np.zeros((len(tenDq), 25))

    for ind in range(len(tenDq)):

        mfc = False
    
        t0 = time.time()
        t00 = t0
        #para, specPara = ReadFromInterface()
        para.printParas()
        N_ligand = para.get_N_ligand()
        nup = para.get_nup()
        ndn = para.get_ndn()
        N = N_ligand + 3 + 5
        
        hop, U_onsite, U_ext, J_ext, U_rest, U_pddp, U_dpdp = SetModelParas(para, tenDq[ind])
        
        HsizeEst = 10000
        Hsp_i = np.zeros((HsizeEst), dtype = np.uint64)
        Hsize_i = GenHsp2p1Hole(N_ligand, nup, ndn, Hsp_i)
        Hsp_i = Hsp_i[:Hsize_i]
        print '    Hilbert Space size with core-hole (intermediate state) is ', Hsp_i.size
        print ' '
        for i in range(Hsp_i.size):
            print i, bin(Hsp_i[i])
        print '    Generating Hilbert space takes ' + str(time.time() - t0) + ' seconds'
        t0 = time.time()
        #n_i = averageElectron(N, Hsp_i, ctempv)    
        sparseH_i = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                    Hsp_i, False, None)


        a = sparseH_i.todense()
        #for i in range(25):
        #    for j in range(25):
        #        if (abs(a[i,j]) > 0.001):
        #            print i, '|', bin(Hsp_0[i])[5:10], bin(Hsp_0[i])[13:18],'>', \
        #                  j, '|', bin(Hsp_0[j])[5:10], bin(Hsp_0[j])[13:18],'>', \
        #                  round(a[i,j].real, 4)


        print ''
        print '    Generating groundstate Hamiltonian takes ' + \
            str(round(time.time() - t0, 3)) + ' seconds'
        t0 = time.time()
            
        E_i, v_i = scipy.sparse.linalg.eigsh(sparseH_i, k=25, M=None, sigma=None, \
            which='SA', v0=None, ncv=40, maxiter=500, tol=1e-07, \
            return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')
        print 'eigenvalues ', E_i


        # test wavefunction

        #for i in range(25):
        #    if(True):
        #        print 'E_0', E_0[i]
        #        print ' v = '
        #        for j in range(25):
        #            if(abs(v_0[j,i]) > 0.001):
        #                print v_0[j,i], '|', bin(Hsp_0[j])[5:10], bin(Hsp_0[j])[13:18],'>', j 
        #        print ' '

        # test wavefunction done

        data[ind, :] = np.sort(E_i[:]) + 3525

    for i in range(25):
        plt.plot(tenDq, data[:,i], 'k') 
    plt.ylim([-7, -3.5])
    plt.show()

    return
 
 
def main():
    para = ModelParas()
    
    # http://www.quanty.org/documentation/tutorials/small_programs_a_quick_start/energy_level_diagram?s[]=tanabe&s[]=sugano
    #=================================
    #--> List of parameters:
    #--> F2dd    = 11.142
    #--> F4dd    =  6.874
    #--> F0dd    = U+(F2dd+F4dd)*2/63 (U = 0.0) = 0.5719
    #--> zeta_3d = 0.081
    #--> Bz      = 0.000001
    #=================================
    #rA = F_0 - 49 * F_4
    #rB = F_2 - 5 * F_4
    #rC = 35 F_4
    
    #F_0 = F^0
    #F_2 = F^2 / 49
    #F_4 = F^4 / 441
    #=================================
    #rC = 35/441 * 6.874 = 0.5455
    #rB = 11.142/49 - 5/441 * 6.874 = 0.14945
    #rA = F0 - 49/441 * 6.874 = (11.142 + 6.874)*2/63 - 49/441 * 6.874 = -0.1918   
    

    #para.setParas(F0 = 0.5719, F2 = 11.142, F4 = 6.874, F0_pd = 0.0, F2_pd = 0.155, \
    #          G1_pd = 0.26688, G3_pd = 0.00928653, N_ligand = 0, E_core = 705.0, \
    #          lambda_SO = 8.5, seedname = 'NiO', \
    #          orders = [5, 6, 1, 7, 4, 3, 2, 8], \
    #          nup = 7, ndn = 7)

    para.setParas(F0 = 0.0, F2 = 0.0, F4 = 0.0, F0_pd = 0.448, F2_pd = 6.667, \
              G1_pd = 4.922, G3_pd = 2.796, N_ligand = 0, E_core = 705.0, \
              lambda_SO = 8.5, seedname = 'NiO', \
              orders = [5, 6, 1, 7, 4, 3, 2, 8], \
              nup = 7, ndn = 7)
              
    test_Tanabe_Sugano(para) 

    return              
 
if __name__ == "__main__":
    main()                 
