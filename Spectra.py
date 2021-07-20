
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

from SetModelParas import SetXASPolarParas, SetModelParas
from GenHsp import GenHspCoreFull, GenHsp2p1Hole
from GenMatrix import GenOStatesSpinUp, GenODaggerStatesSpinUpDown, GenMatrix, averageElectron
from GenMatrix import GenODaggerStatesSpinUp, GenODaggerStatesSpinDown
from CmplxContFracExpan import CmplxContFracExpan
from BiCGS import BiCGS
#from ReadFromInterface import ReadFromInterface
#from ClassDefine import ModelParas, SpectraParas

#para, specPara = ReadFromInterface()

def XAS(para, specPara):

    mfc = False

    t0 = time.time()
    t00 = t0
    #para, specPara = ReadFromInterface()
    para.printParas()
    N_ligand = para.get_N_ligand()
    nup = para.get_nup()
    ndn = para.get_ndn()
    N = N_ligand + 3 + 5
    
    hop, U_onsite, U_ext, J_ext, U_rest, U_pddp, U_dpdp = SetModelParas(para)
    
    HsizeEst = 10000000
    Hsp_0 = np.zeros((HsizeEst), dtype = np.uint64)
    Hsize_0 = GenHspCoreFull(N_ligand, nup, ndn, Hsp_0)
    Hsp_0 = Hsp_0[:Hsize_0]
    #print 'Hsp[0]', bin(Hsp_0[0])
    #for i in range(Hsp_0.size):
    #    print i, bin(Hsp_0[i])
    print '    Generating Hilbert space takes ' + str(time.time() - t0) + ' seconds'
    print '    Hilbert Space size w/o core-hole (ground state) is ', Hsp_0.size
    print ' '
    
    
    sparseH_0 = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                          Hsp_0, False, None)
    
    print ''
    print '    Generating groundstate Hamiltonian takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()
        
    E_0, v_0 = scipy.sparse.linalg.eigsh(sparseH_0, k=2, M=None, sigma=None, \
        which='SA', v0=None, ncv=10, maxiter=500, tol=1e-07, \
        return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')
    v_0 = v_0[:,0]
    E_0 = E_0[0]
    print ''
    print '    Groundstate energy', E_0, '(eV)'    
    print '    Diagonalizing Hamiltonian matrices takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()



    n_0 = averageElectron(N, Hsp_0, v_0)    
    print n_0, sum(n_0)
    if(mfc):
        # Mean-Field correction    
        print n_0, sum(n_0)   
        sparseH_0 = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                              Hsp_0, True, n_0)
        E_0, v_0 = scipy.sparse.linalg.eigsh(sparseH_0, k=2, M=None, sigma=None, \
            which='SA', v0=None, ncv=10, maxiter=500, tol=1e-07, \
            return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')
        v_0 = v_0[:,0]
        E_0 = E_0[0]
        print ''
        print '    New groundstate energy', E_0, '(eV)'    
        print '    Mean-field correction takes ' + \
            str(round(time.time() - t0, 3)) + ' seconds'
        t0 = time.time()
        n_0p = averageElectron(N, Hsp_0, v_0)    
        print n_0p, sum(n_0p)   
        # Mean-Field correction done
                                                                                        
 
                      
    Hsp_i = np.zeros((HsizeEst), dtype = np.uint64)
    Hsize_i = GenHsp2p1Hole(N_ligand, nup, ndn, Hsp_i)
    Hsp_i = Hsp_i[:Hsize_i]
    print '    Hilbert Space size with core-hole (intermediate state) is ', Hsp_i.size
    print ' '
    #for i in range(Hsp_i.size):
    #    print i, bin(Hsp_i[i])
    print '    Generating Hilbert space takes ' + str(time.time() - t0) + ' seconds'
    t0 = time.time()
    #n_i = averageElectron(N, Hsp_i, ctempv)    
    sparseH_i = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                Hsp_i, mfc, n_0)    

    print ''
    print '    Generating Intermediatestate Hamiltonian takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()
    
    # Mean-Field correction for Intermediate State
    E_i, v_i = scipy.sparse.linalg.eigsh(sparseH_i, k=3, M=None, sigma=None, \
        which='LM', v0=None, ncv=10, maxiter=400, tol=1e-07, \
        return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')
    v_i = v_i[:,0]
    E_i = E_i[0]
    print ''
    #for i in range(10):
    #    print 'Intermediate state eigenvector ', v_i[i]

    print '    Intermediate state minimum energy', E_i, '(eV)'   
    n_i = averageElectron(N, Hsp_i, v_i)    
    print n_i, sum(n_i)    
    '''
    sparseH_i = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                          Hsp_i, True, n_i)
    print '' 
    print '    Mean-field correction takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()
    '''
    
            
    specX = np.linspace(specPara.start_Ein, specPara.end_Ein, \
                        num = specPara.div_Ein + 1, endpoint = True)
    specY = np.zeros((specX.size))                    
                                                                                                                                                                            
    for pol in range(len(specPara.theta_in)):                                   

        tsxas = time.time()
        Gn = SetXASPolarParas(specPara.theta_in[pol], specPara.phi_in[pol])
        print 'Gn ', Gn
        print ''
        print '    SetXASPolarParas takes ' + \
            str(round(time.time() - tsxas, 3)) + ' seconds'

        tgos = time.time()
        ctempv = GenOStatesSpinUp(v_0, Hsp_0, Hsp_i, Gn, N)
        print '    GenOSttesSpinUp takes ' + \
            str(round(time.time() - tgos, 3)) + ' seconds'

        tcnorm = time.time()
        print 'norm ctempv ', np.vdot(ctempv, ctempv)
        n_i = averageElectron(N, Hsp_i, ctempv)  
        print 'n_i and sum ', n_i, sum(n_i)
        print '    averageElectron takes ' + \
            str(round(time.time() - tcnorm, 3)) + ' seconds'
        '''
        sparseH_i = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                Hsp_i, True, n_i)    
        print '' 
        print '    Mean-field correction takes ' + \
            str(round(time.time() - t0, 3)) + ' seconds'
        '''   
        tcfe = time.time()             
        myspecY = CmplxContFracExpan(ctempv, E_0, sparseH_i, specX, specPara.epsilone_Ein)    
        specY += myspecY 
        print ''
        print '    Continued fraction expansion takes ' + \
            str(round(time.time() - tcfe, 3)) + ' seconds for pol no.', pol
        t0 = time.time()
    
    filename = os.getcwd() + '/' + para.get_seedname() + '.XAS.dat'
    #f = open(filename, "w")
    #for mm in range(specY.size):
    #    f.write(str(round(specX[mm],3)) + ' ' + str(round(specY[mm],6)) + '\n')
    specXsize = specX.size
    specX = specX.reshape((specXsize,1))
    specY = specY.reshape((specXsize,1))
    np.savetxt(filename, np.concatenate((specX,specY),axis=1))
    #f.close()
        
    print ''
    print '    Calculated XAS for' + para.get_seedname() + \
        'has been written at file: ' + para.get_seedname() + '.XAS.dat'    
    print ''    
    print '    Finishing the program takes ' + \
        str(round(time.time() - t00, 3)) + ' seconds'

    return
 
       

def RIXS(para, specPara):
    
    mfc = False
    t0 = time.time()
    t00 = t0

    # --------------------------------------------------------------------------
        
    # Step 1 Generate Model parameters 
        
    #para, specPara = ReadFromInterface()
    para.printParas()
    N_ligand = para.get_N_ligand()
    nup = para.get_nup()
    ndn = para.get_ndn()
    N = N_ligand + 3 + 5
    
    hop, U_onsite, U_ext, J_ext, U_rest, U_pddp, U_dpdp = SetModelParas(para)

    
         
    # --------------------------------------------------------------------------
         
    # Step 2 Generating Hilbert Space         
                
    HsizeEst = comb(2*N-6, nup+ndn-6, exact=True)
    #HsizeEst = 10000000
    Hsp_0 = np.zeros((HsizeEst), dtype = np.uint64)
    Hsize_0 = GenHspCoreFull(N_ligand, nup, ndn, Hsp_0)
    Hsp_0 = Hsp_0[:Hsize_0]
    print '    Hilbert Space size w/o core-hole (ground state) is ', Hsp_0.size
    print ' '


    #N_core = 3
    #N_d = 5
    #N = N_d + N_ligand + N_core
    HsizeEst = 6 * comb(2*N-6, nup+ndn-5, exact=True)
    Hsp_i = np.zeros((HsizeEst), dtype = np.uint64)
    Hsize_i = GenHsp2p1Hole(N_ligand, nup, ndn, Hsp_i)
    Hsp_i = Hsp_i[:Hsize_i]
    print '    Hilbert Space size with core-hole (intermediate state) is ', Hsp_i.size
    print ' '
    
    # non spin flip
    HsizeEst = comb(2*N-6, nup+ndn-6, exact=True)
    Hsp_f_spin0 = np.zeros((HsizeEst), dtype = np.uint64)
    Hsize_f_spin0 = GenHspCoreFull(N_ligand, nup, ndn, Hsp_f_spin0)
    Hsp_f_spin0 = Hsp_f_spin0[:Hsize_f_spin0]

    Hsp_f_spin2 = np.zeros((HsizeEst), dtype = np.uint64)
    Hsize_f_spin2 = GenHspCoreFull(N_ligand, nup+1, ndn-1, Hsp_f_spin2)
    Hsp_f_spin2 = Hsp_f_spin2[:Hsize_f_spin2]
    print '    Hilbert Space size w/o core-hole (final state) non-spin-flip is ', Hsp_f_spin0.size
    print '    Hilbert Space size w/o core-hole (final state) for spin-flip is ', Hsp_f_spin2.size
    print ' '
    print '    Generating Hilbert space takes ' + str(time.time() - t0) + ' seconds'
    t0 = time.time()
         
         
         
    # --------------------------------------------------------------------------
         
    # Step 3 Generating Hamiltonian Matrix 
         
         
    sparseH_0 = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                          Hsp_0, False, None)
    print ''
    print '    Generating groundstate Hamiltonian matrices takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()

    E_0, v_0 = scipy.sparse.linalg.eigsh(sparseH_0, k=2, M=None, sigma=None, \
        which='SA', v0=None, ncv=5, maxiter=200, tol=1e-05, \
        return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')
    v_0 = v_0[:,0]
    E_0 = E_0[0]
    print ''
    print '    Groundstate energy', E_0, '(eV)'
    print '    Diagonalizing Hamiltonian matrices takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()



    n_0 = averageElectron(N, Hsp_0, v_0)
    # Mean-Field correction    -------------------------------------------------
    if(mfc): 
        print n_0, sum(n_0)
        sparseH_0 = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                              Hsp_0, True, n_0)
        E_0, v_0 = scipy.sparse.linalg.eigsh(sparseH_0, k=2, M=None, sigma=None, \
            which='SA', v0=None, ncv=5, maxiter=200, tol=1e-05, \
            return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')
        v_0 = v_0[:,0]
        E_0 = E_0[0]
        print ''
        print '    New groundstate energy', E_0, '(eV)'
        print '    Mean-field correction takes ' + \
            str(round(time.time() - t0, 3)) + ' seconds'
        t0 = time.time()
        n_0p = averageElectron(N, Hsp_0, v_0)
        print n_0p, sum(n_0p)
        print 'sparseH_0 col', sparseH_0.tocoo().col
    
    # Mean-Field correction done   ---------------------------------------------
    
                    
            
    sparseH_i = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                          Hsp_i, mfc, n_0)
    print ''
    print '    Generating intermediatestate Hamiltonian matrices takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()
    
    
    sparseH_f_spin0 = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                          Hsp_f_spin0, mfc, n_0)
    sparseH_f_spin2 = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                          Hsp_f_spin2, mfc, n_0)
    print ''
    print '    Generating finalstate Hamiltonian matrices takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()
            
            
            

    # --------------------------------------------------------------------------
         
    # Step 4 Calculate RIXS cross-section
                
    filename = os.getcwd() + '/' + para.get_seedname() + '.RIXS.dat'
    #f = open(filename, "w")

        
    '''
    ! y-pol
    theta_in = 0.0
    phi_in = pi/2
    ! x-pol
    theta_out = 0.0
    phi_out = 0.0
    !do phi_out = 0, pi, pi/10
        !write(*,*) 'theta, phi at ', theta, phi
    '''
                        
    for eneY in range(int(specPara.div_Ein) + 1): 
        
        specX = np.linspace(specPara.start_Eloss, specPara.end_Eloss, \
                        num = specPara.div_Eloss + 1, endpoint = True)
        specY_spin0 = np.zeros((int(specPara.div_Eloss) + 1))        
        specY_spin2 = np.zeros((int(specPara.div_Eloss) + 1))        
        #specY[:] = 0.0
        

        Ein = eneY*1.0/(specPara.div_Ein)*(specPara.end_Ein-specPara.start_Ein) + \
            specPara.start_Ein
        z = Ein + E_0 + specPara.epsilone_Ein * 1j
        
        t0 = time.time()   
        
        for pol_in in range(len(specPara.theta_in)):
            for pol_out in range(len(specPara.theta_out)):
                    
                Gn = SetXASPolarParas(specPara.theta_in[pol_in], specPara.phi_in[pol_in])
                ctempv = GenOStatesSpinUp(v_0, Hsp_0, Hsp_i, Gn, N)
            
                # generate sparseH_i_z matrix (sparseH_i - z)
                
                scoo = sparseH_i.tocoo()
                Hvalue_z = np.zeros((Hsize_i), dtype = np.complex128)
                IndexI_z = np.zeros((Hsize_i), dtype = np.intc)
                IndexJ_z = np.zeros((Hsize_i), dtype = np.intc)
                for i in range(Hsize_i):
                    IndexI_z[i] = i
                    IndexJ_z[i] = i
                    Hvalue_z[i] = -1 * z
                sparseH_i_z = scipy.sparse.coo_matrix((np.append(scoo.data, Hvalue_z),\
                        (np.append(scoo.col, IndexI_z), np.append(scoo.row, IndexJ_z))),\
                        shape=(Hsize_i,Hsize_i)).tocsr()
              
                ctempv, info = scipy.sparse.linalg.bicgstab(sparseH_i_z, ctempv, x0=None, \
                    tol=1e-06, maxiter=50, xtype=None, M=None, callback=None)
                #call BiCGS(z, Hsize_i, SprSize_i, ctempv, IndexI_i, IndexJ_i, sparseH_i)
            
                #print 'ctempv ', ctempv

                Gn = SetXASPolarParas(specPara.theta_out[pol_out], specPara.phi_out[pol_out])

                pvector_spin0 = GenODaggerStatesSpinUpDown(ctempv, Hsp_i, Hsp_f_spin0, Gn, N)
                myspecY = CmplxContFracExpan(pvector_spin0, E_0, sparseH_f_spin0, \
                                             specX, specPara.epsilone_Eloss)
                specY_spin0 += myspecY

                pvector_spin2 = GenODaggerStatesSpinUpDown(ctempv, Hsp_i, Hsp_f_spin2, Gn, N)
                myspecY = CmplxContFracExpan(pvector_spin2, E_0, sparseH_f_spin2, \
                                             specX, specPara.epsilone_Eloss)
                specY_spin2 += myspecY
                
        print '        Incident Energy ', Ein, 'eV, BiCGSTAB finishes within ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
            
        #for mm in range(specY.size):
        #    f.write(Ein, specX[mm], specY[mm])
        specXsize = specX.size
        specX = specX.reshape((specXsize,1))
        specY_spin0 = specY_spin0.reshape((specXsize,1))
        specY_spin2 = specY_spin2.reshape((specXsize,1))
        Einarray = np.zeros((specXsize,1))
        Einarray[:,:] = Ein

        f_handle_spin0 = file(os.getcwd() + '/' + para.get_seedname() + '.RIXS-nonspinflip.dat', 'a')
        np.savetxt(f_handle_spin0, np.concatenate((Einarray, specX, specY_spin0),axis=1))
        f_handle_spin0.close()
            
        f_handle_spin2 = file(os.getcwd() + '/' + para.get_seedname() + '.RIXS-spinflip.dat', 'a')
        np.savetxt(f_handle_spin2, np.concatenate((Einarray, specX, specY_spin2),axis=1))
        f_handle_spin2.close()

    #f.close()
    
    print ''
    print '    Calculated RIXS for' + para.get_seedname() + \
        'has been written at file:' + filename  
    print ''            
    print '    Finishing the program takes ' + str(round(time.time() - t00, 3)) \
        + ' seconds'
    return


def RIXS2(para, specPara):
   
    mfc = False
    t0 = time.time()
    t00 = t0

    # --------------------------------------------------------------------------
    
    # Step 1 Generate Model parameters 
    
    #para, specPara = ReadFromInterface()
    para.printParas()
    N_ligand = para.get_N_ligand()
    nup = para.get_nup()
    ndn = para.get_ndn()
    N = N_ligand + 3 + 5
   
    hop, U_onsite, U_ext, J_ext, U_rest, U_pddp, U_dpdp = SetModelParas(para)
   
    
    # --------------------------------------------------------------------------
    
    # Step 2 Generating Hilbert Space         
    
    HsizeEst = comb(2*N-6, nup+ndn-6, exact=True)
    #HsizeEst = 10000000
    Hsp_0 = np.zeros((HsizeEst), dtype = np.uint64)
    Hsize_0 = GenHspCoreFull(N_ligand, nup, ndn, Hsp_0)
    Hsp_0 = Hsp_0[:Hsize_0]
    #for i in range(Hsp_0.size):
    #    print i, bin(Hsp_0[i])
    print '    Generating Hilbert space takes ' + str(time.time() - t0) + ' seconds'
    print '    Hilbert Space size w/o core-hole (ground state) is ', Hsp_0.size
    print ' '


    #N_core = 3
    #N_d = 5
    #N = N_d + N_ligand + N_core
    HsizeEst = 6 * comb(2*N-6, nup+ndn-5, exact=True)
    Hsp_i = np.zeros((HsizeEst), dtype = np.uint64)
    Hsize_i = GenHsp2p1Hole(N_ligand, nup, ndn, Hsp_i)
    Hsp_i = Hsp_i[:Hsize_i]
    print '    Hilbert Space size with core-hole (intermediate state) is ', Hsp_i.size
    print ' '
    #for i in range(Hsp_i.size):
    #    print i, bin(Hsp_i[i])

    # non spin flip
    HsizeEst = comb(2*N-6, nup+ndn-6, exact=True)
    Hsp_f_spin0 = np.zeros((HsizeEst), dtype = np.uint64)
    Hsize_f_spin0 = GenHspCoreFull(N_ligand, nup, ndn, Hsp_f_spin0)
    Hsp_f_spin0 = Hsp_f_spin0[:Hsize_f_spin0]

    Hsp_f_spin2 = np.zeros((HsizeEst), dtype = np.uint64)
    Hsize_f_spin2 = GenHspCoreFull(N_ligand, nup+1, ndn-1, Hsp_f_spin2)
    Hsp_f_spin2 = Hsp_f_spin2[:Hsize_f_spin2]
    print '    Hilbert Space size w/o core-hole (final state) non-spin-flip is ', Hsp_f_spin0.size
    print '    Hilbert Space size w/o core-hole (final state) for spin-flip is ', Hsp_f_spin2.size
    print ' '
    print '    Generating Hilbert space takes ' + str(time.time() - t0) + ' seconds'
    t0 = time.time()



    # --------------------------------------------------------------------------

    # Step 3 Generating Hamiltonian Matrix 

    sparseH_0 = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                          Hsp_0, False, None)
    print ''
    print '    Generating groundstate Hamiltonian matrices takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()

    E_0, v_0 = scipy.sparse.linalg.eigsh(sparseH_0, k=20, M=None, sigma=None, \
        which='SA', v0=None, ncv=30, maxiter=600, tol=1e-08, \
        return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')

    # analyze average electron for each excited state    
   # for i in range(20):
   #     f_handle_n = file(os.getcwd() + '/' + para.get_seedname() + '.averagen-nsf.dat', 'a')
   #     n_0 = averageElectron(N, Hsp_0, v_0[:,i])
   #     top = np.append(E_0[i] - E_0[0], n_0)
   #     np.savetxt(f_handle_n, top.reshape(1, 23))
   #     f_handle_n.close()



    v_0 = v_0[:,0]
    E_0 = E_0[0]
    print ''
    print '    Groundstate energy', E_0, '(eV)'
    print '    Diagonalizing Hamiltonian matrices takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()



    n_0 = averageElectron(N, Hsp_0, v_0)
    print n_0, sum(n_0)
    # Mean-Field correction    -------------------------------------------------
    if(mfc):
        print n_0, sum(n_0)
        sparseH_0 = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                              Hsp_0, True, n_0)
        E_0, v_0 = scipy.sparse.linalg.eigsh(sparseH_0, k=2, M=None, sigma=None, \
            which='SA', v0=None, ncv=10, maxiter=600, tol=1e-08, \
            return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')
        v_0 = v_0[:,0]
        E_0 = E_0[0]
        print ''
        print '    New groundstate energy', E_0, '(eV)'
        print '    Mean-field correction takes ' + \
            str(round(time.time() - t0, 3)) + ' seconds'
        t0 = time.time()
        n_0p = averageElectron(N, Hsp_0, v_0)
        print n_0p, sum(n_0p)
        print 'sparseH_0 col', sparseH_0.tocoo().col

    # Mean-Field correction done   ---------------------------------------------



    sparseH_i = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                          Hsp_i, mfc, n_0)
    print ''
    print '    Generating intermediatestate Hamiltonian matrices takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()


    sparseH_f_spin0 = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                          Hsp_f_spin0, mfc, n_0)
    sparseH_f_spin2 = GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
                          Hsp_f_spin2, mfc, n_0)
    print ''
    print '    Generating finalstate Hamiltonian matrices takes ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
    t0 = time.time()




    # --------------------------------------------------------------------------

    # Step 4 Calculate RIXS cross-section

    #f = open(filename, "w")


    '''
    ! y-pol
    theta_in = 0.0
    phi_in = pi/2
    ! x-pol
    theta_out = 0.0
    phi_out = 0.0
    !do phi_out = 0, pi, pi/10
        !write(*,*) 'theta, phi at ', theta, phi
    '''

    f_handle_spin0 = file(os.getcwd() + '/' + para.get_seedname() + '.RIXS-nsf.dat', 'w')
    f_handle_spin0.close()
    f_handle_spin2 = file(os.getcwd() + '/' + para.get_seedname() + '.RIXS-sf.dat', 'w')
    f_handle_spin2.close()


    for eneY in range(int(specPara.div_Ein) + 1):

        this_epsilone_Ein = 0.05+float(eneY)/specPara.div_Ein * specPara.epsilone_Ein
        #this_epsilone_Ein = specPara.epsilone_Ein
        
        specX = np.linspace(specPara.start_Eloss, specPara.end_Eloss, \
                        num = specPara.div_Eloss + 1, endpoint = True)
        specY_spin0 = np.zeros((int(specPara.div_Eloss) + 1))
        specY_spin2 = np.zeros((int(specPara.div_Eloss) + 1))
        #specY[:] = 0.0


        Ein = eneY*1.0/(specPara.div_Ein)*(specPara.end_Ein-specPara.start_Ein) + \
            specPara.start_Ein
        z = Ein + E_0 + specPara.epsilone_Ein * 1j
        #z = Ein + E_0 + this_epsilone_Ein * 1j
        #print 'eneY = ', eneY, ' Ein = ', Ein, 'eV, this_epsilone_Ein = ', this_epsilone_Ein, ' eV'
        #print 'specPara.div_Ein = ', specPara.div_Ein, 'specPara.epsilone_Ein = ', specPara.epsilone_Ein

        # generate sparseH_i_z matrix (sparseH_i - z)
        #scoo = sparseH_i.tocoo()
        #Hvalue_z = np.zeros((Hsize_i), dtype = np.complex64)
        #IndexI_z = np.zeros((Hsize_i), dtype = np.intc)
        #IndexJ_z = np.zeros((Hsize_i), dtype = np.intc)
        #for i in range(Hsize_i):
        #    IndexI_z[i] = i
        #    IndexJ_z[i] = i
        #    Hvalue_z[i] = -1 * z
        #sparseH_i_z = scipy.sparse.coo_matrix((np.append(scoo.data, Hvalue_z),\
        #        (np.append(scoo.col, IndexI_z), np.append(scoo.row, IndexJ_z))),\
        #        shape=(Hsize_i,Hsize_i)).tocsr()

        t0 = time.time()

        for pol_in in range(len(specPara.theta_in)):
            for pol_out in range(len(specPara.theta_out)):

                Gn = SetXASPolarParas(specPara.theta_in[pol_in], specPara.phi_in[pol_in])
                ctempv = GenOStatesSpinUp(v_0, Hsp_0, Hsp_i, Gn, N)

                #ctempvp, info = scipy.sparse.linalg.bicgstab(sparseH_i_z, ctempv, x0=None, \
                #    tol=1e-05, maxiter=50, xtype=None, M=None, callback=None)
                ctempvp, info, tsum = BiCGS(z, ctempv, sparseH_i, 200, 0.00001)

                if(info == -1):
                    print '            bicgstab fails at incoming energy ', Ein, ' tsum = ', tsum

                Gn = SetXASPolarParas(specPara.theta_out[pol_out], specPara.phi_out[pol_out])


                pvector = GenODaggerStatesSpinUp(ctempvp, Hsp_i, Hsp_0, Gn, N)
                #myspecY = CmplxContFracExpan(pvector_spin0, E_0, sparseH_0, \
                #                             specX, specPara.epsilone_Eloss, 100)
                myspecY = np.zeros((int(specPara.div_Eloss) + 1))
                for eneX in range(specX.size):

                    this_epsilone_Eloss = 0.05+float(eneX)/specX.size * specPara.epsilone_Eloss
                    #zp = E_0 + specPara.epsilone_Eloss * 1j + specX[eneX]
                    zp = E_0 + this_epsilone_Eloss * 1j + specX[eneX]

                    ctempvpp, info, tsum = BiCGS(zp, pvector, sparseH_f_spin0, 200, 0.00001)
                    #if(info == -1):
                    #    print '            bicgstab fails at incoming ', Ein, \
                    #          ' and energy loss ', specX[eneX], ' tsum = ', tsum
                    myspecY[eneX] += (np.vdot(ctempvpp, pvector)).imag / 3.1416
                specY_spin0 += myspecY
                #print 'myspecY', myspecY



                pvector_spin2 = GenODaggerStatesSpinDown(ctempvp, Hsp_i, Hsp_f_spin2, Gn, N)
                #myspecY = CmplxContFracExpan(pvector_spin2, E_0, sparseH_f_spin2, \
                #                             specX, specPara.epsilone_Eloss)
                myspecY1 = np.zeros((int(specPara.div_Eloss) + 1))
                for eneX in range(specX.size):

                    this_epsilone_Eloss = 0.05+float(eneX)/specX.size * specPara.epsilone_Eloss
                    zp = E_0 + this_epsilone_Eloss * 1j + specX[eneX]

                    ctempvpp, info, tsum = BiCGS(zp, pvector_spin2, sparseH_f_spin2, 200, 0.00001)
                    #if(info == -1):
                    #    print '            bicgstab fails at incoming ', Ein, \
                    #          ' and energy loss ', specX[eneX], ' tsum = ', tsum
                    myspecY1[eneX] += (np.vdot(ctempvpp, pvector_spin2)).imag / 3.1416
                    #print '         energy loss ', specX[eneX], 'eV with broadening', this_epsilone_Eloss
                specY_spin2 += myspecY1

        print '        Incident Energy ', Ein, 'eV, BiCGSTAB finishes within ' + \
        str(round(time.time() - t0, 3)) + ' seconds'
        print '           broadening ', specPara.epsilone_Ein

        #for mm in range(specY.size):
        #    f.write(Ein, specX[mm], specY[mm])
        specXsize = specX.size
        specX = specX.reshape((specXsize,1))
        specY_spin0 = -1 * specY_spin0.reshape((specXsize,1))
        specY_spin2 = -1 * specY_spin2.reshape((specXsize,1))
        Einarray = np.zeros((specXsize,1))
        Einarray[:,:] = Ein

        filename = os.getcwd() + '/' + para.get_seedname() + '.RIXS-nsf.dat'
       
        f_handle_spin0 = file(os.getcwd() + '/' + para.get_seedname() + '.RIXS-nsf.dat', 'a')
        np.savetxt(f_handle_spin0, np.concatenate((Einarray, specX, specY_spin0),axis=1))
        f_handle_spin0.close()

        f_handle_spin2 = file(os.getcwd() + '/' + para.get_seedname() + '.RIXS-sf.dat', 'a')
        np.savetxt(f_handle_spin2, np.concatenate((Einarray, specX, specY_spin2),axis=1))
        f_handle_spin2.close()

    #f.close()

    print ''
    print '    Calculated RIXS for' + para.get_seedname() + \
        'has been written at file:' + filename
    print ''
    print '    Finishing the program takes ' + str(round(time.time() - t00, 3)) \
        + ' seconds'
    return

