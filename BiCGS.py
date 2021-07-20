'''
June 17, 2018 Change from fortran code to python  


'''
import numpy as np
import time

def BiCGS(z, v, sparseH, niter_CG, tol_CG):


    ##def BiCGS(z, Hsizet, SprSizet, groundH, IndexI, IndexJ, sparseH)
    #use NumOfOrbitalAndElectrons; use ConstantParas; use ScanRegion
    #implicit none
    #
    #double complex, INTENT(IN) :: z
    #integer*8, INTENT(IN) :: SprSizet, Hsizet
    #double complex, DIMENSION(Hsizet), INTENT(INOUT) :: groundH
    #integer*8, DIMENSION(SprSizet), INTENT(IN) :: IndexI, IndexJ
    #double complex, DIMENSION(SprSizet), INTENT(IN) :: sparseH
    #
    #double complex, allocatable:: x0(:), xn(:)
    #double complex, allocatable:: r0(:), rn(:), r00(:)
    #double complex, allocatable:: p0(:), pn(:)
    #double complex, allocatable:: v0(:), vn(:)
    #double complex, allocatable:: s(:), at(:)
    #double complex, allocatable:: dvtemp(:)
    #double complex :: alpha, rho0, rhon
    #double complex :: myrhon, mytsum, tsum, mytsum1, tsum1
    #double complex :: betha, omega0, omegan
    #integer*8 :: spr_size, nloc, ii, jj
    #double precision :: time1, time2
    #     !z = (0.0d0, 0.0d0)
    #   !****************        STEP 3.A        *****************
    #   !********              Initialization           **********
    #
    #nloc = Hsizet
    #spr_size = SprSizet
    #
    #allocate(x0(nloc),xn(nloc))
    #allocate(p0(nloc),pn(nloc))
    #allocate(v0(nloc),vn(nloc))
    #allocate(r0(nloc),rn(nloc),r00(nloc))
    #allocate(s(nloc),at(nloc))
    #allocate(dvtemp(nloc))

    Hsize = v.size
    x0 = np.zeros((Hsize), dtype=np.complex128)
    xn = np.zeros((Hsize), dtype=np.complex128)
    p0 = np.zeros((Hsize), dtype=np.complex128)
    pn = np.zeros((Hsize), dtype=np.complex128)
    v0 = np.zeros((Hsize), dtype=np.complex128)
    vn = np.zeros((Hsize), dtype=np.complex128)
    r0 = np.zeros((Hsize), dtype=np.complex128)
    rn = np.zeros((Hsize), dtype=np.complex128)
    r00 = np.zeros((Hsize), dtype=np.complex128)
    s = np.zeros((Hsize), dtype=np.complex128)
    at = np.zeros((Hsize), dtype=np.complex128)
    dvtemp = np.zeros((Hsize), dtype=np.complex128)

    r0[:] = v[:]
    r00[:] = r0[:]
    # Question? For such non-array variables, how to set double complex?
    rho0 = 1.0 + 0.0j
    alpha = 1.0 + 0.0j
    omega0 = 1.0 + 0.0j
    tsum = 0.0

    t0 = time.time()
    for ii in range(niter_CG):  

        rhon = np.vdot(r00,r0)
        betha = (rhon/rho0)*(alpha/omega0)
        pn = r0 + betha*(p0-omega0*v0)
        vn = sparseH.dot(pn) - z * pn  # Not sure if the order is correct, we will see   
        alpha = rhon / np.vdot(r00,vn)
        s = r0 - alpha*vn
        at = sparseH.dot(s) - z * s
        omegan = np.vdot(at,s) / np.vdot(at,at)
        xn = x0 + alpha * pn + omegan * s
        rn = s - omegan * at
        dvtemp = sparseH.dot(xn) - z * xn - v
        tsum = np.vdot(dvtemp, dvtemp)
        #print '      inside CG ', ii, tsum
        if(tsum.real < tol_CG):
            t2 = time.time()
            #print '      CG finishes at:', t2 - t0, 'secs for precision ', tsum
            return xn, ii, tsum
     
        x0[:] = xn[:] 
        rho0 = rhon
        p0[:] = pn[:]
        r0[:] = rn[:]
        omega0 = omegan
        v0[:] = vn[:]

    return xn, -1, tsum

