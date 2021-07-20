'''
!!Cheng-Chien Chen

!!Date:02/27/2011.
!!Lanczos code serial version. Note: careful about renormalizaton of Lanczos coefficients.

!!Date:05/30/2011.
!!Modifying the Lacnzos code to perform continued fraction expansion (serial).
!!The serial version is capable of dealing with Hsize=10^7 and obtain spectra in hours.

!!Note:the input Hamiltonian matirx (sparseH.txt) and OGS.txt are now complex.
'''
#from scipy import sparse
import numpy as np
import time

def CmplxContFracExpan(v, E_0, sparseH, specX, epsilone_CFE):

    t0 = time.time()
    iter_num = 100
    #print 'iter_num ', iter_num
    Hsize = v.size
    #print 'v.size and Hsize ', v.size
    #epsilone_CFE = 0.3
    specY = np.zeros((specX.size))
    
    carray_1 = np.zeros((Hsize), dtype = np.complex128)
    phi_n0 = np.zeros((Hsize), dtype = np.complex128)
    phi_n1 = np.zeros((Hsize), dtype = np.complex128)
    phi_nm = np.zeros((Hsize), dtype = np.complex128)
    
    #phi_n0[:Hsize] = v[:Hsize]
    phi_n0[:] = v[:]
    
    work_an = np.zeros((iter_num))
    work_bn = np.zeros((iter_num))


#!!Specify the order of continued fraction expansions.
#!!Note: Usually iter_num=50 already gives truncation error < 10-6;
#!!      For iter_num > 100 orthogonality may leak and result in "NaN".
#!!************************************************
    
    for iiter in range(iter_num):
    
        #print 'iiter ', iiter
        #!!Constructing H|phi_n0>: later want to do in one step 
        # so only 3 vectors are needed.
        carray_1 = sparseH.dot(phi_n0)
        #print 'phi_n0 size ', phi_n0.size
        #print 'carray_1 size ', carray_1.size
        
        #Calculate <phi_n0|H|phi_n0>:
        ctemp1= np.vdot(phi_n0, carray_1)
        rtemp1= np.vdot(phi_n0, phi_n0)
        a_n=ctemp1/rtemp1
        
        if (iiter == 0): b_n=0.0
        phi_n1[:] = carray_1[:] - a_n*phi_n0[:] - b_n*phi_nm[:]
        
        #print 'a_n, b_n ', a_n, b_n
        work_an[iiter]=np.real(a_n)
        work_bn[iiter]=np.real(b_n)
        
        rtemp2 = np.vdot(phi_n1, phi_n1)
        b_n=rtemp2/rtemp1
        
        phi_nm[:] = phi_n0[:]
        phi_n0[:] = phi_n1[:]
        phi_n1[:] = 0.0
        
        phi_n0[:] = phi_n0[:]/np.sqrt(rtemp2)
        phi_nm[:] = phi_nm[:]/np.sqrt(rtemp2)
    
    carray_1[:Hsize] = v[:Hsize]
    rtemp1 = np.vdot(carray_1, carray_1)
    
    t1 = time.time()
    #print ' Part 1 of CFE takes ', str(round(t1 - t0, 3)), 'seconds'


    for Eind in range(specX.size):
        
        ctemp1 =  specX[Eind] + E_0 + epsilone_CFE * 1j
        ctemp4 = 0.0
        alist = range(1,iter_num)
        #do jj=iter_num, 2, -1
        for jj in alist[::-1]:
            ctemp3 = work_bn[jj]/(ctemp1-work_an[jj] - ctemp4)
            ctemp4 = ctemp3
            
        con_ct1 = rtemp1/(ctemp1-work_an[0]-ctemp4)
        specY[Eind] = np.imag(con_ct1)/(-1*np.pi)

    #call CPU_time(time2)
    #print 'Complex CFE (secs):', time2-time1

    t2 = time.time()
    #print ' Part 2 of CFE takes ', str(round(t2 - t1, 3)), 'seconds'

    return specY

