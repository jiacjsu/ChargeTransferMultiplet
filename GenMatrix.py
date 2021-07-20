
# format of Gn:
#        write(*,*) 'Polarization Gn'
#         write(*,*) Gn(1,3), Gn(2,3), Gn(3,3)
#         write(*,*) Gn(4,3), Gn(5,3)
# H_e complex array:

import numpy as np
import scipy
import scipy.sparse
import ctypes

def sumEveryBit(tempa):
    #cdef int tempsum=0
    tempsum = 0
    tempa1 = tempa
    while tempa1:
        tempsum += tempa1 % 2
        tempa1 = tempa1 >> 1
    return tempsum
 

def GenOStatesSpinUp(H_0, Hsp_src, Hsp_dst, Gn, N):
    
    N_core = 3
    N_3d = 5
    H_e = np.zeros((Hsp_dst.size), dtype = np.complex128)

    lib = ctypes.cdll['./GenMatrixCore.so']
    GenOStatesCoreSpinUp = lib['GenOStatesCoreSpinUp']

    Gn_f = Gn.flatten('C')

    GenOStatesCoreSpinUp(ctypes.c_void_p(H_0.ctypes.data), len(Hsp_src), ctypes.c_void_p(Hsp_src.ctypes.data), \
                   ctypes.c_void_p(H_e.ctypes.data), len(Hsp_dst), ctypes.c_void_p(Hsp_dst.ctypes.data), \
                   ctypes.c_void_p(Gn_f.ctypes.data), N)

    #print ' test H_e sum inside python subroutin ', sum(H_e[:])
    #print ' test H_e norm inside python ........ ', np.vdot(H_e, H_e)
    return H_e

   
#def GenOStates(H_0, Hsp_src, Hsp_dst, Gn, N):
#    
#    N_core = 3
#    N_3d = 5
#    H_e = np.zeros((Hsp_dst.size), dtype = np.complex64)
#
#    #cdef int jj, popo, toto
#    for jj in range(Hsp_src.size):
#        statup = Hsp_src[jj] / (2**N)
#        statdn = Hsp_src[jj] % (2**N)
#        for popo in range(N-N_core, N):
#            #if (.not.BTEST(statdn,popo)): cycle   !c2p, dn
#            for toto in range(N_3d):
#                maskpo = 1 << popo
#                maskto = 1 << toto
#                if (statup & maskpo): 
#                    istatup = statup & (~maskpo)
#                    tempsign = (-1)**(sumEveryBit(statup >> (popo+1)))
#                    if not(istatup & maskto):
#                        fstatup = istatup | maskto
#                        tempsign *= (-1)**(sumEveryBit(istatup >> (toto+1)))                     
#                               
#                        #lBinarySearch(Hsp_dst,Hsize_dst,statup,fstatdn,l)
#                        l = np.searchsorted(Hsp_dst, fstatup * (2**N) + statdn)
#                        
#                        H_e[l] += Gn[toto, popo - (N-N_core)]*tempsign*H_0[jj]
#
#    return H_e


def GenODaggerStatesSpinUpDown(H_0, Hsp_src, Hsp_dst, Gn, N):
    

    N_core = 3
    N_3d = 5
    H_e = np.zeros((Hsp_dst.size), dtype = np.complex128)

    lib = ctypes.cdll['./GenMatrixCore.so']
    GenODaggerStatesCoreSpinUpDown = lib['GenODaggerStatesCoreSpinUpDown']

    Gn_f = Gn.flatten('C')

    GenODaggerStatesCoreSpinUpDown(ctypes.c_void_p(H_0.ctypes.data), len(Hsp_src), ctypes.c_void_p(Hsp_src.ctypes.data), \
                   ctypes.c_void_p(H_e.ctypes.data), len(Hsp_dst), ctypes.c_void_p(Hsp_dst.ctypes.data), \
                   ctypes.c_void_p(Gn_f.ctypes.data), N)

    #print ' test H_e sum inside python subroutin ', sum(H_e[:])
    #print ' test H_e norm inside python ........ ', np.vdot(H_e, H_e)
    return H_e



def GenODaggerStatesSpinUp(H_0, Hsp_src, Hsp_dst, Gn, N):

    N_core = 3
    N_3d = 5
    H_e = np.zeros((Hsp_dst.size), dtype = np.complex128)

    lib = ctypes.cdll['./GenMatrixCore.so']
    GenODaggerStatesCoreSpinUp = lib['GenODaggerStatesCoreSpinUp']

    Gn_f = Gn.flatten('C')

    GenODaggerStatesCoreSpinUp(ctypes.c_void_p(H_0.ctypes.data), len(Hsp_src), ctypes.c_void_p(Hsp_src.ctypes.data), \
                   ctypes.c_void_p(H_e.ctypes.data), len(Hsp_dst), ctypes.c_void_p(Hsp_dst.ctypes.data), \
                   ctypes.c_void_p(Gn_f.ctypes.data), N)

    #print ' test H_e sum inside python subroutin ', sum(H_e[:])
    #print ' test H_e norm inside python ........ ', np.vdot(H_e, H_e)
    return H_e


def GenODaggerStatesSpinDown(H_0, Hsp_src, Hsp_dst, Gn, N):

    N_core = 3
    N_3d = 5
    H_e = np.zeros((Hsp_dst.size), dtype = np.complex128)

    lib = ctypes.cdll['./GenMatrixCore.so']
    GenODaggerStatesCoreSpinDown = lib['GenODaggerStatesCoreSpinDown']

    Gn_f = Gn.flatten('C')

    GenODaggerStatesCoreSpinDown(ctypes.c_void_p(H_0.ctypes.data), len(Hsp_src), ctypes.c_void_p(Hsp_src.ctypes.data), \
                   ctypes.c_void_p(H_e.ctypes.data), len(Hsp_dst), ctypes.c_void_p(Hsp_dst.ctypes.data), \
                   ctypes.c_void_p(Gn_f.ctypes.data), N)

    #print ' test H_e sum inside python subroutin ', sum(H_e[:])
    #print ' test H_e norm inside python ........ ', np.vdot(H_e, H_e)
    return H_e





def averageElectron(N, Hsp, v):
    n = np.zeros((2*N))
    Hsize = Hsp.size
    lib = ctypes.cdll['./GenMatrixCore.so']
    averageElectronCore = lib['averageElectronCore']
    averageElectronCore(N, Hsize, ctypes.c_void_p(Hsp.ctypes.data), \
                        ctypes.c_void_p(v.ctypes.data),  \
                        ctypes.c_void_p(n.ctypes.data))
    return n




#def averageElectron(N, Hsp, v):
#    n = np.zeros((2*N))
#    Hsize = Hsp.size
#    N_eg = 2
#    N_t2g = 3
#    N_d = 5
#    for i in range(Hsize):
#        statupdn = Hsp[i]
#        for j in range(2*N):
#            if statupdn & (1<<j):
#                n[j] += np.real(np.vdot(v[i], v[i]))
#    for j in range(N):
#        t = (n[j] + n[j+N])/2.0
#        n[j] = t
#        n[j+N] = t
#    tsum = sum(n[:N_eg]) / N_eg
#    for i in range(N_eg) + range(N, N+N_eg):
#        n[i] = tsum      
#    tsum = sum(n[N_eg:N_d]) / N_t2g
#    for i in range(N_eg, N_d) + range(N + N_eg, N+N_d):
#        n[i] = tsum  
#    modv = np.real(np.vdot(v, v))
#    n = n/modv
#    return n
    

# corr is the tag for Mean field correction
# v_0 is the ground state wavefunction
def GenMatrix(N, hop, U_onsite, U_ext, J_ext, U_pddp, U_dpdp, U_rest, \
              Hsp, corr, n_0): 
    print ''
    print ''
    #print 'U_c ', U_c
    Hsize = Hsp.size
    
    E_corr = MFCorr(U_onsite, U_ext, N, corr, n_0)
    if corr: 
        print 'Mean-Field Correction'
        print E_corr[:N]
        print E_corr[N:]
    
    nsparse = Hsize*60
    print('   Hsize', Hsize)
    print('   Allocate IndexI, Hvalue size', nsparse)
  #  for i in range(Hsize):
  #      print('   ', "{0:b}".format(Hsp[i]))
  #  for i in range(N):
  #      for j in range(N):
  #          print('    hop',i,j, hop[i,j] )
  #  for i in range(N):
  #      for j in range(N):
  #          print('    J_ext',i,j, J_ext[i,j] )
  #  for i in range(N):
  #      for j in range(N):
  #          print('    U_ext',i,j, U_ext[i,j] )
  #  for i in range(5):
  #      for j in range(3):
  #          for k in range(5):
  #              for l in range(3):
  #                  print('    U_dpdp',i,j,k,l, U_dpdp[i,j,k,l] )


    Hvalue = np.zeros((nsparse), dtype = np.complex128)
    IndexI = np.zeros((nsparse), dtype = np.intc)
    IndexJ = np.zeros((nsparse), dtype = np.intc)


    lib = ctypes.cdll['./GenMatrixCore.so']
    GenMatrixCore = lib['GenMatrixCore']

    hop_f = hop.flatten('C')
    U_ext_f = U_ext.flatten('C')
    J_ext_f = J_ext.flatten('C')
    U_rest_f = U_rest.flatten('C')
    U_pddp_f = U_pddp.flatten('C')
    U_dpdp_f = U_dpdp.flatten('C')


    nsparse = GenMatrixCore(Hsize, N, \
                            ctypes.c_void_p(hop_f.ctypes.data), \
                            ctypes.c_void_p(U_onsite.ctypes.data), \
                            ctypes.c_void_p(U_ext_f.ctypes.data), \
                            ctypes.c_void_p(J_ext_f.ctypes.data), \
                            ctypes.c_void_p(U_rest_f.ctypes.data), \
                            ctypes.c_void_p(U_pddp_f.ctypes.data), \
                            ctypes.c_void_p(U_dpdp_f.ctypes.data), \
                            ctypes.c_void_p(Hsp.ctypes.data), \
                            corr, ctypes.c_void_p(E_corr.ctypes.data), \
                            ctypes.c_void_p(Hvalue.ctypes.data), \
                            ctypes.c_void_p(IndexI.ctypes.data), \
                            ctypes.c_void_p(IndexJ.ctypes.data))

    Hvalue = Hvalue[:nsparse]
    IndexI = IndexI[:nsparse]
    IndexJ = IndexJ[:nsparse]

    #for i in range(len(Hvalue)):
    #    print 'Matrix ', IndexI[i], IndexJ[i], Hvalue[i]

   # C conterpart
   # int GenMatrixCore(int Hsize, int N, double complex *hop, double *U_onsite,
   #                 double *U_ext, double *J_ext, double *U_pddp, double *U_dpdp,
   #                 unsigned int *Hsp, bool corr, double *E_corr,
   #                 double complex *Hvalue, int *IndexI, int *IndexJ){


   # for i in range(Hsize):
   #     statupdn = Hsp[i]
   #     tdict = {}
   #     d = HopHami(statupdn, hop, N, corr, E_corr)
   #     for k in d:
   #         if k in tdict:
   #             tdict[k] += d[k]
   #         else:
   #             tdict[k] = d[k]
   #     d = UonsiteHami(statupdn, U_onsite, N)
   #     for k in d:
   #         if k in tdict:
   #             tdict[k] += d[k]
   #         else:
   #             tdict[k] = d[k]
   #     d = UextHami(statupdn, U_ext, N)
   #     for k in d:
   #         if k in tdict:
   #             tdict[k] += d[k]
   #         else:
   #             tdict[k] = d[k]      
   #     d = JextHami(statupdn, J_ext, N)
   #     for k in d:
   #         if k in tdict:
   #             tdict[k] += d[k]
   #         else:
   #             tdict[k] = d[k]       
   #     d = dpMulti(statupdn, U_pddp, U_dpdp, N)
   #     for k in d:
   #         if k in tdict:
   #             tdict[k] += d[k]
   #         else:
   #             tdict[k] = d[k]       
   #    # d = ddMulti(statupdn, U_rest, N)
   #    # for k in d:
   #    #     if k in tdict:
   #    #         tdict[k] += d[k]
   #    #     else:
   #    #         tdict[k] = d[k] 
   #             
   #     # combine the same index
   #     #Hvalue = np.array(Hvalue)
   #     #IndexI = np.array([i] * len(Hindex))
   #     #IndexJ = np.array(Hindex)
   #     Hvalue = []
   #     Hindex = []
   #     for k in tdict:
   #         Hvalue.append(tdict[k])
   #         Hindex.append(np.searchsorted(Hsp, k))
   #     Hvalue = np.array(Hvalue)
   #     IndexI = np.array([i] * len(Hindex))
   #     IndexJ = np.array(Hindex)
   #     if i%1000 == 0:
   #         #print '        tdict ', tdict
   #         print '        i, Hsize, min(IndexJ), max(IndexJ), size of dict ', \
   #         i, Hsize, min(IndexJ), max(IndexJ), len(tdict)
    sparseH = scipy.sparse.coo_matrix((Hvalue,(IndexI,IndexJ)),\
              shape=(Hsize,Hsize)).tocsr()
   #     if i == 0:
   #         sparseH = sparseH_t
   #     else:
   #         sparseH += sparseH_t
        
    return sparseH
    


def MFCorr(U_onsite, U_ext, N, corr, n):

    E_corr = np.zeros((2*N))

    #cdef int pp,kk,j    
    #cdef float Uijcorrect
    if corr:
        for pp in range(2*N):
            kk = (pp+N) % (2*N)
            Uijcorrect = 0.0
            for j in range(2*N):
                if (j % N) != (pp % N):
                    Uijcorrect += U_ext[pp % N, j % N] * n[j]
            E_corr[pp] = -0.5 * U_onsite[pp%N] * n[kk] - 0.5 * Uijcorrect
    t2g = 0.0
    N_d = 5
    for pp in range(2,N_d) + range(N+2, N+N_d):
        t2g += E_corr[pp]
    t2g = t2g/(2*3)
    for pp in range(2,N_d) + range(N+2, N+N_d):
        E_corr[pp] = t2g

    eg = 0.0
    for pp in range(2) + range(N, N+2):
        eg += E_corr[pp]
    eg = eg/(2*2)
    for pp in range(2) + range(N, N+2):
        E_corr[pp] = eg
    return E_corr


 
# Refer to Ilkyu's notes on MF correction    
def MFCorr_d(U_onsite, U_ext, N, corr, n):
   
    E_corr = np.zeros((2*N))
    
    #cdef int pp,kk,j    
    #cdef float Uijcorrect
    if corr:
        for pp in range(2*N):
            kk = (pp+N) % (2*N)
            Uijcorrect = 0.0
            for j in range(2*N):
                if (j % N) != (pp % N):
                    Uijcorrect += U_ext[pp % N, j % N] * n[j]
            E_corr[pp] = -0.5 * U_onsite[pp%N] * n[kk] - 0.5 * Uijcorrect
    t_d = 0.0
    N_d = 5
    for pp in range(N_d):
        t_d += (E_corr[pp] + E_corr[pp+N])
    t_d = t_d / 10.0
    for pp in range(N_d):
        E_corr[pp] = t_d 
        E_corr[pp+N] = t_d 
    #for pp in range(2,N_d) + range(N+2, N+N_d):
    #    t2g += E_corr[pp]
    #t2g = t2g/(2*3)
    #for pp in range(2,N_d) + range(N+2, N+N_d): 
    #    E_corr[pp] = t2g   
    #                                  
    #eg = 0.0
    #for pp in range(2) + range(N, N+2):
    #    eg += E_corr[pp]
    #eg = eg/(2*2)
    #for pp in range(2) + range(N, N+2): 
    #    E_corr[pp] = eg
    return E_corr     
          
