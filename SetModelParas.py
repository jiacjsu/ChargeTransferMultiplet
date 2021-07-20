
# Chunjing Jia Oct 02, 2017
#    calculate the paramters for the transition metal 3d or 4d systems with 2p core 

# problem 1: for calculating 2-particle interactions, U and J, what are the index range?
# problem 2: output parameters, how to pass?
# problem 3: input parameters calculated from Python script


# input parameters:
#
#   rA, rB, rB:         Racah paramters for 3d electrons
#   F_0, F_2, G_1, G_3: Condon-Shortley parameters (originally defined as the uppercase,
#                       but here used as the lowercase) for 2p 3d electrons
#      			core-hole potential U_2p3d is embeded in F_0
#   lambda_SO:          the spin orbital coupling constant for 2p core levels
#   seedname:		the filename for the _hr.dat hopping integrals
#   orders:		since the wannier downfolding has different orders than defined 
#   			in this file, we need to change th eordering

# output parameters:
#
#   hop[2*ntol, 2*ntol]:	all on-site energies and hopping integrals read from _hr.dat file
#   U_d:	 		diagonal terms of Coulomb interactions for 3d electrons
#   U_ext[ntol, ntol]:	 	off diagonal terms of Coulomb interactions for 3d electrons
#   J_ext[ntol, ntol]:		exchange terms for 3d electrons
#   U_rest[5,5,5,5]: 			other multiplet interaction terms for 3d electrons
#   U_pddp[3,5,5,3], U_dpdp[5,3,5,3]: 		multiplet terms for 2p 3d electrons 

import numpy as np
import os
#import ClassDefine 

def SetXASPolarParas(theta, phi):
    #'''
    #!   ! z-pol
    #!   theta_in = pi/2
    #!   phi_in = 0.0
    #!   ! y-pol
    #!   theta_out = 0.0
    #!   phi_out = pi/2
    #!   ! x-pol
    #!   theta_out = 0.0
    #!   phi_out = 0.0
    #'''
    alpha_p = np.cos(theta)*np.cos(phi)
    beta_p  = np.cos(theta)*np.sin(phi)
    gamma_p = -1*np.sin(theta)
    Gn = np.zeros((5,3))
    Gn[0,0] = -1*alpha_p
    Gn[0,1] = -1*beta_p
    Gn[0,2] = 2.0*gamma_p;
    Gn[1,0] = np.sqrt(3.0)*alpha_p
    Gn[1,1] = -np.sqrt(3.0)*beta_p
    Gn[2,0] = np.sqrt(3.0)*beta_p
    Gn[2,1] = np.sqrt(3.0)*alpha_p
    Gn[3,0] = np.sqrt(3.0)*gamma_p
    Gn[3,2] = np.sqrt(3.0)*alpha_p
    Gn[4,1] = np.sqrt(3.0)*gamma_p
    Gn[4,2] = np.sqrt(3.0)*beta_p
    
    #Gn[0,0] = -0.060515*alpha_p
    #Gn[0,1] = -0.060515*beta_p
    #Gn[0,2] = 0.60515*gamma_p
    #Gn[1,0] = 0.15722*alpha_p
    #Gn[1,1] = -0.15722*beta_p
    #Gn[2,0] = 0.15722*beta_p
    #Gn[2,1] = 0.15722*alpha_p
    #Gn[3,0] = 0.20963*gamma_p
    #Gn[3,2] = 0.20963*alpha_p
    #Gn[4,1] = 0.20963*gamma_p
    #Gn[4,2] = 0.20963*beta_p
    
    #for i in range(5):
    #    for j in range(3):
    #        if abs(Gn[i,j]) < 0.02:
    #            Gn[i,j] = 0.02
    #Gn[:,:] = 1.0
    return Gn
    

def SetModelParas(para):

    rA = para.get_rA()
    rB = para.get_rB()
    rC = para.get_rC()
    F_0 = para.get_F_0()
    F_2 = para.get_F_2()
    G_1 = para.get_G_1()
    G_3 = para.get_G_3()
    N_ligand = para.get_N_ligand()
    E_core = para.get_E_core()
    lambda_SO = para.get_lambda_SO()
    seedname = para.get_seedname()
    orders = para.get_orders()
    #para.printParas()
    
    N_d = 5
    N_core = 3
    ntol = N_d + N_ligand + N_core
    ntolmax = 28
    #print 'ntol ', ntol

    # index order: N_d, N_ligand, L-edge
    # orders in this Python codes 
    # d orbitals: 1=d_{3z^2-r^2}; 2=d_{x^2-y^2}; 3: dxy; 4: dxz; 5: dyz.
    # p orbitals: 1 px, 2 py, 3 pz 
    # need to transfer between Wannier output and this convention

    #'''
    #small index represents LSB
    #index 0 through ntol-1: spin down
    #index ntol through 2*ntol-1: spin up
    #'''
    hop = np.zeros((2*ntolmax, 2*ntolmax), dtype=np.complex128)
    tdim = 0
    f = open(os.getcwd() + '/' + seedname + '_hr.dat', "r")
    f.readline()
    f.readline()
    line = f.readline()
    #print int(line.split()[0])/15
    for i in range(int(line.split()[0])/15 + 1):
        f.readline()
    while True:
        line = f.readline()
        #print line
        if line:
            nx, ny, nz, oa, ob, val_real, val_imag = line.split()[:]
            oa = int(oa)
            ob = int(ob)
            val_real = float(val_real)
            val_imag = float(val_imag)
            if np.sqrt(val_real**2 + val_imag**2) > abs(hop[orders[oa-1]-1, orders[ob-1]-1]):
                hop[orders[oa-1]-1, orders[ob-1]-1] = val_real + val_imag * 1j
        else:
            break
    for i in range(ntol - N_core, ntol):
        hop[i,i] -= E_core

    #TenDq = 1.050
    #Hop[0,0] = 0.6*tenDq
    #Hop[1,1] = 0.6*tenDq 

    #Hop[2,2] = -0.4*tenDq
    #Hop[3,3] = -0.4*tenDq
    #Hop[4,4] = -0.4*tenDq
    #Eshiftpi = -200.0
    #for i in range(5, 5+6):
    #    hop[i,i] -= 10
    #for i in range(5+6, 5+7):
    #    hop[i,i] += 10
    
    # electron-hole language change
    #for i in range(ntol):
    #    hop[i,i] = -1 * hop[i,i]

    #Eshiftsigma = -40.0
    #for i in range(8, 10):
    #    hop[i,i] += Eshiftsigma
    
    #print 'Delta shift for pi ligand is ', Eshiftpi, ' eV'    
    #print 'Delta shift for sigma ligand is ', Eshiftsigma, ' eV'    
    print 'On-site Energy: ', 'dz2', hop[0,0], 'dx2-y2', hop[1,1]
    print '                ', 'dxy', hop[2,2], 'dxz', hop[3,3], 'dyz', hop[4,4]
    print '                ', 'ligands ', hop[5,5], hop[6,6], hop[7,7], hop[8,8], hop[9,9], hop[10,10], hop[11,11]
    print '                ', 'cores', hop[12,12], hop[13,13], hop[14,14]
    #print '                ', 'cores', hop[8,8], hop[9,9], hop[10,10]
     
    #for i in range(ntol):
    #    for j in range(ntol):
    #        if(hop[i,j] != 0.0):
    #            print ' print hop i, j, ', i, j, hop[i,j]
    print ' '
    #print ' '
   
    '''
    ! Defining the 2p S-O matrix:
    ! U_SO_2p(3 by 2 by 3 by 2).
    ! Spin-orbit coupling strength = xi_2p*U_SO_2p
    ! the second and fourth index: 2 stands up, 1 stands down
    ! the first and third index: 1 px, 2 py, 3 pz
    '''
    U_SO_2p = np.zeros((3,2,3,2), dtype=np.complex128)
    U_SO_2p[0,0,1,0]= 0.5j
    U_SO_2p[0,0,2,1]= -0.5
    U_SO_2p[0,1,1,1]= -0.5j
    U_SO_2p[0,1,2,0]= 0.5
    U_SO_2p[1,0,0,0]= -0.5j 
    U_SO_2p[1,0,2,1]= -0.5j
    U_SO_2p[1,1,0,1]= 0.5j   
    U_SO_2p[1,1,2,0]= -0.5j
    U_SO_2p[2,0,0,1]= 0.5
    U_SO_2p[2,0,1,1]= 0.5j
    U_SO_2p[2,1,0,0]= -0.5
    U_SO_2p[2,1,1,0]= 0.5j
    for i in range(ntol):
        for j in range(ntol):
            hop[i+ntol, j+ntol] = hop[i,j]
    for i in range(3):
        for j in range(3):
            hop[i+N_d+N_ligand,j+N_d+N_ligand] += U_SO_2p[i,0,j,0] * lambda_SO
            hop[i+N_d+N_ligand,j+N_d+N_ligand+ntol] += U_SO_2p[i,0,j,1] * lambda_SO
            hop[i+N_d+N_ligand+ntol,j+N_d+N_ligand] += U_SO_2p[i,1,j,0] * lambda_SO
            hop[i+N_d+N_ligand+ntol,j+N_d+N_ligand+ntol] += U_SO_2p[i,1,j,1] * lambda_SO

    '''
    !*****Coulomb Interaction Parameters*****
    
    !Hund's couplings (Coulomb Interaction) for the 5 d-orbitals
    !1. The onsite repulstion U, which has been calculated before.
    !2. Coulomb repulsion U_ext for different orbitas, all spins (diagonal elements)
    !3. The interband exchange J, for all spins, different orbitals
    !4. J_ext for different spins (pair hopping), different orbitals
    '''

    U_d = rA+4.0*rB+3.0*rC 
    U_onsite = np.zeros((ntolmax))
    U_onsite[:5] = U_d
    print 'U_onsite ', U_d
    #U = np.zeros((2*ntol, 2*ntol))
    U_ext = np.zeros((ntolmax, ntolmax))
    J_ext = np.zeros((ntolmax, ntolmax), dtype = np.float_)

    U_ext[0,1]= rA-4.00*rB+rC    # f0 - 4f2 + 6f4,  2112
    U_ext[0,2]= rA-4.00*rB+rC    #                  2332
    U_ext[1,0]= rA-4.00*rB+rC    # 
    U_ext[2,0]= rA-4.00*rB+rC

    U_ext[0,3]= rA+2.00*rB+rC    # f0 + 2f2 - 24f4, 2442
    U_ext[0,4]= rA+2.00*rB+rC    #                  2552
    U_ext[3,0]= rA+2.00*rB+rC
    U_ext[4,0]= rA+2.00*rB+rC

    U_ext[1,2]= rA+4.00*rB+rC    # f0 + 4f2 - 34f4, 1331
    U_ext[2,1]= rA+4.00*rB+rC

    U_ext[3,2]= rA-2.00*rB+rC    # f0 - 2f2 - 4f4, 3443
    U_ext[2,3]= rA-2.00*rB+rC
    U_ext[2,4]= rA-2.00*rB+rC    #                 3553  
    U_ext[4,2]= rA-2.00*rB+rC

    U_ext[3,4]= rA-2.00*rB+rC    # f0 - 2f2 - 4f4, 4554
    U_ext[4,3]= rA-2.00*rB+rC

    U_ext[1,3]= rA-2.00*rB+rC    # f0 - 2f2 - 4f4, 1441
    U_ext[1,4]= rA-2.00*rB+rC    #                 1551
    U_ext[3,1]= rA-2.00*rB+rC
    U_ext[4,1]= rA-2.00*rB+rC

    J_ext[0,1]= 4.00*rB+rC       # 4f2 + 15f4, 2121
    J_ext[1,0]= 4.00*rB+rC
    J_ext[0,2]= 4.00*rB+rC       #           , 2323
    J_ext[2,0]= 4.00*rB+rC

    J_ext[0,3]= rB+rC            # f2 + 30f4, 2424
    J_ext[3,0]= rB+rC           
    J_ext[0,4]= rB+rC            #          , 2525
    J_ext[4,0]= rB+rC

    J_ext[1,2]= rC               # 35f4, 1313
    J_ext[2,1]= rC

    J_ext[1,3]= 3.00*rB+rC       # 3f2 + 20f4, 1414
    J_ext[3,1]= 3.00*rB+rC
    J_ext[1,4]= 3.00*rB+rC       #           , 1515
    J_ext[4,1]= 3.00*rB+rC

    J_ext[2,3]= 3.00*rB+rC       #           , 3434
    J_ext[2,4]= 3.00*rB+rC       #           , 3535
    J_ext[3,2]= 3.00*rB+rC
    J_ext[3,4]= 3.00*rB+rC       #           , 4545
    J_ext[4,2]= 3.00*rB+rC
    J_ext[4,3]= 3.00*rB+rC

    #U_ext = np.zeros((ntolmax, ntolmax))
    #J_ext = np.zeros((ntolmax, ntolmax), dtype = np.float_)

#    for i in range(5):
#        for j in range(5):
#            U[i,j] = U_ext[i,j] / 2.0
#            U[i,j+ntol] = U_ext[i,j] / 2.0
#            U[i+ntol,j] = U_ext[i,j] / 2.0
#            U[i+ntol,j+ntol] = U_ext[i,j] / 2.0
#    for i in range(5):
#	U[i,i+ntol] = U_d / 2.0
#	U[i+ntol,i] = U_d / 2.0
#
    '''
    do ii=1,5
    do jj=1,5
          U(ii,jj) = U_ext(ii,jj)/2
          U(ii+N,jj) = U_ext(ii,jj)/2
          U(ii,jj+N) = U_ext(ii,jj)/2
          U(ii+N,jj+N) = U_ext(ii,jj)/2
    enddo
    enddo
    do jj=1,5
          U(jj+N,jj) = U_d/2
          U(jj,jj+N) = U_d/2
    enddo
    '''
    '''
    !********* Exchange Interaction  *********
    '''
    '''
    do ii=1, 5
       do jj=1, 5
          if(U_ext(ii,jj).ne.U_ext(jj,ii)) then
             write(*,*) 'non-hermitian U_ext:',ii,jj
             stop
          endif
    
          if(J_ext(ii,jj).ne.J_ext(jj,ii)) then
             write(*,*) 'non-hermitian J_ext:',ii,jj
             stop
          endif
       enddo
    enddo
    '''
    '''
    !*********Four particle terms U_rest, U_pddp, U_dpdp **********
    U_multiplet = 0.0
    U_rest=0.0d0 U_pddp=0.0d0; U_dpdp=0.0d0;
    
    !!Note:
    !!F_2=F^2/35 G_1=G^1/15; G_3=G^3/245;
    
    !Ti(4+): Sawatzky
    !F_0=0.0d0 F_2=0.18d0; G_1=0.308d0; G_3=2.63d0/245;
    !F_2=F_2*0.8d0
    !G_1=G_1*0.8d0
    !G_3=G_3*0.8d0
    
    
    !!!Fe(3+): Sawatzky et al., PRB 42, 5459)
    !!!Fe(2+): Sawatzky et al., PRB 42, 5459)
    !!F_0=0.0d0 F_2=0.155; G_1=0.26688; G_3=0.00928653;
    '''
    
    #U_multiplet = np.zeros((ntol, ntol, ntol, ntol))
    U_rest = np.zeros((N_d, N_d, N_d, N_d))
    
    U_rest[1,0,3,3]=  np.sqrt(3.00)*rB   # sqrt(3)f2 - 5*sqrt(3)f4, 4412
    U_rest[1,3,0,3]= np.sqrt(3.00)*rB
    U_rest[3,0,3,1]=  np.sqrt(3.00)*rB
    U_rest[3,3,0,1]=  np.sqrt(3.00)*rB
    U_rest[0,3,1,3]= np.sqrt(3.00)*rB
    U_rest[0,1,3,3]=  np.sqrt(3.00)*rB
    U_rest[3,3,1,0]=  np.sqrt(3.00)*rB
    U_rest[3,1,3,0]= np.sqrt(3.00)*rB

    U_rest[1,0,4,4]= -np.sqrt(3.00)*rB   # -sqrt(3)f2 + 5*sqrt(3)f4, 5512
    U_rest[1,4,0,4]= -np.sqrt(3.00)*rB
    U_rest[4,0,4,1]= -np.sqrt(3.00)*rB
    U_rest[4,4,0,1]= -np.sqrt(3.00)*rB
    U_rest[0,4,1,4]=  -np.sqrt(3.00)*rB
    U_rest[0,1,4,4]= -np.sqrt(3.00)*rB
    U_rest[4,1,4,0]= -np.sqrt(3.00)*rB
    U_rest[4,4,1,0]= -np.sqrt(3.00)*rB

    U_rest[2,0,3,4]=  np.sqrt(3.00)*rB   # sqrt(3)f2 - 5*sqrt(3)f4, 2354
    U_rest[2,3,0,4]= np.sqrt(3.00)*rB
    U_rest[4,0,3,2]=  np.sqrt(3.00)*rB
    U_rest[4,3,0,2]=  np.sqrt(3.00)*rB
    U_rest[0,2,4,3]=  np.sqrt(3.00)*rB
    U_rest[0,4,2,3]=  np.sqrt(3.00)*rB
    U_rest[3,2,4,0]= np.sqrt(3.00)*rB
    U_rest[3,4,2,0]=  np.sqrt(3.00)*rB

    U_rest[2,0,4,3]=  np.sqrt(3.00)*rB  # 2345
    U_rest[2,4,0,3]= np.sqrt(3.00)*rB
    U_rest[3,0,4,2]=  np.sqrt(3.00)*rB
    U_rest[3,4,0,2]=  np.sqrt(3.00)*rB
    U_rest[0,2,3,4]=  np.sqrt(3.00)*rB
    U_rest[0,3,2,4]= np.sqrt(3.00)*rB
    U_rest[4,2,3,0]= np.sqrt(3.00)*rB
    U_rest[4,3,2,0]=  np.sqrt(3.00)*rB

    U_rest[4,0,2,3]=  -2.00*np.sqrt(3.00)*rB
    U_rest[4,2,0,3]= -2.00*np.sqrt(3.00)*rB
    U_rest[3,0,2,4]=  -2.00*np.sqrt(3.00)*rB
    U_rest[3,2,0,4]= -2.00*np.sqrt(3.00)*rB
    U_rest[0,4,3,2]=  -2.00*np.sqrt(3.00)*rB  # -2*sqrt(3)f2 + 10*sqrt(3)f4, 2453
    U_rest[0,3,4,2]= -2.00*np.sqrt(3.00)*rB
    U_rest[2,3,4,0]= -2.00*np.sqrt(3.00)*rB
    U_rest[2,4,3,0]= -2.00*np.sqrt(3.00)*rB

    U_rest[4,0,1,4]=   2.00*np.sqrt(3.00)*rB # 2*sqrt(3)f2 - 10*sqrt(3)f4, 5215
    U_rest[4,1,0,4]=  2.00*np.sqrt(3.00)*rB
    U_rest[1,4,4,0]= 2.00*np.sqrt(3.00)*rB
    U_rest[0,4,4,1]=  2.00*np.sqrt(3.00)*rB

    U_rest[3,1,0,3]= -2.00*np.sqrt(3.00)*rB  # -2*sqrt(3)f2 + 10*sqrt(3)f4, 4214
    U_rest[3,0,1,3]=  -2.00*np.sqrt(3.00)*rB  #added
    U_rest[0,3,3,1]= -2.00*np.sqrt(3.00)*rB
    U_rest[1,3,3,0]= -2.00*np.sqrt(3.00)*rB

    U_rest[2,1,3,4]=  3.00*rB  # 3*f2 - 15*f4, 1354
    U_rest[2,3,1,4]= 3.00*rB
    U_rest[4,1,3,2]=  3.00*rB
    U_rest[4,3,1,2]=  3.00*rB
    U_rest[1,2,4,3]=  3.00*rB
    U_rest[1,4,2,3]= 3.00*rB
    U_rest[3,2,4,1]= 3.00*rB
    U_rest[3,4,2,1]=  3.00*rB

    U_rest[2,1,4,3]= -3.00*rB  # -3*f2 + 15*f4, 1345
    U_rest[2,4,1,3]= -3.00*rB
    U_rest[3,1,4,2]= -3.00*rB
    U_rest[3,4,1,2]=  -3.00*rB
    U_rest[1,2,3,4]= -3.00*rB
    U_rest[1,3,2,4]= -3.00*rB
    U_rest[4,2,3,1]= -3.00*rB
    U_rest[4,3,2,1]=  -3.00*rB

    #U_rest = np.zeros((N_d, N_d, N_d, N_d))

    U_pddp = np.zeros((3,5,5,3))
    U_dpdp = np.zeros((5,3,5,3))

    U_pddp[0,0,0,0]= -(-F_0+2.00*F_2)
    U_pddp[1,0,0,1]= -(-F_0+2.00*F_2)
    U_pddp[2,0,0,2]= -(-F_0-4.00*F_2)
    U_pddp[0,1,1,0]= -(-F_0-2.00*F_2)
    U_pddp[1,1,1,1]= -(-F_0-2.00*F_2)
    U_pddp[2,1,1,2]= -(-F_0+4.00*F_2)
    U_pddp[0,2,2,0]= -(-F_0-2.00*F_2)
    U_pddp[1,2,2,1]= -(-F_0-2.00*F_2)
    U_pddp[2,2,2,2]= -(-F_0+4.00*F_2)
    U_pddp[0,3,3,0]= -(-F_0-2.00*F_2)
    U_pddp[1,3,3,1]= -(-F_0+4.00*F_2)
    U_pddp[2,3,3,2]= -(-F_0-2.00*F_2)
    U_pddp[0,4,4,0]= -(-F_0+4.00*F_2)
    U_pddp[1,4,4,1]= -(-F_0-2.00*F_2)
    U_pddp[2,4,4,2]= -(-F_0-2.00*F_2)

    U_pddp[0,0,1,0]= -(2.00*np.sqrt(3.00)*F_2)
    U_pddp[0,0,2,1]= -(2.00*np.sqrt(3.00)*F_2)
    U_pddp[0,0,3,2]= -(-np.sqrt(3.00)*F_2)
    U_pddp[1,0,1,1]= -(-2.00*np.sqrt(3.00)*F_2)
    U_pddp[1,0,2,0]= -(2.00*np.sqrt(3.00)*F_2)
    U_pddp[1,0,4,2]= -(-np.sqrt(3.00)*F_2)
    U_pddp[2,0,3,0]= -(-np.sqrt(3.00)*F_2)
    U_pddp[2,0,4,1]= -(-np.sqrt(3.00)*F_2)
    U_pddp[0,1,0,0]= -(2.00*np.sqrt(3.00)*F_2)
    U_pddp[0,1,3,2]= -(-3.00*F_2)
    U_pddp[1,1,0,1]= -(-2.00*np.sqrt(3.00)*F_2)
    U_pddp[1,1,4,2]= -(3.00*F_2)
    U_pddp[2,1,3,0]= -(-3.00*F_2)
    U_pddp[2,1,4,1]= -(3.00*F_2)
    U_pddp[0,2,0,1]= -(2.00*np.sqrt(3.00)*F_2)
    U_pddp[0,2,4,2]= -(-3.00*F_2)
    U_pddp[1,2,0,0]= -(2.00*np.sqrt(3.00)*F_2)
    U_pddp[1,2,3,2]= -(-3.00*F_2)
    U_pddp[2,2,3,1]= -(-3.00*F_2)
    U_pddp[2,2,4,0]= -(-3.00*F_2)
    U_pddp[0,3,0,2]= -(-np.sqrt(3.00)*F_2)
    U_pddp[0,3,1,2]= -(-3.00*F_2)
    U_pddp[0,3,4,1]= -(-3.00*F_2)
    U_pddp[1,3,2,2]= -(-3.00*F_2)
    U_pddp[1,3,4,0]= -(-3.00*F_2)
    U_pddp[2,3,0,0]= -(-np.sqrt(3.00)*F_2)
    U_pddp[2,3,1,0]= -(-3.00*F_2)
    U_pddp[2,3,2,1]= -(-3.00*F_2)
    U_pddp[0,4,2,2]= -(-3.00*F_2)
    U_pddp[0,4,3,1]= -(-3.00*F_2)
    U_pddp[1,4,0,2]= -(-np.sqrt(3.00)*F_2)
    U_pddp[1,4,1,2]= -(3.00*F_2)
    U_pddp[1,4,3,0]= -(-3.00*F_2)
    U_pddp[2,4,0,1]= -(-np.sqrt(3.00)*F_2)
    U_pddp[2,4,1,1]= -(3.00*F_2)
    U_pddp[2,4,2,0]= -(-3.00*F_2)

    U_dpdp[0,0,0,0]= -(-G_1-18.00*G_3)
    U_dpdp[0,0,1,0]= -(np.sqrt(3.00)*G_1+3.00*np.sqrt(3.00)*G_3)
    U_dpdp[0,0,2,1]= -(np.sqrt(3.00)*G_1+3.00*np.sqrt(3.00)*G_3)
    U_dpdp[0,0,3,2]= -(-2.00*np.sqrt(3.00)*G_1+9.00*np.sqrt(3.00)*G_3)
    U_dpdp[1,0,0,0]= -(np.sqrt(3.00)*G_1+3.00*np.sqrt(3.00)*G_3)
    U_dpdp[1,0,1,0]= -(-3.00*G_1-24.00*G_3)
    U_dpdp[1,0,2,1]= -(3.00*G_1-21.00*G_3)
    U_dpdp[1,0,3,2]= -(-15.00*G_3)
    U_dpdp[2,0,0,1]= -(np.sqrt(3.00)*G_1+3.00*np.sqrt(3.00)*G_3)
    U_dpdp[2,0,1,1]= -(-3.00*G_1+21.00*G_3)
    U_dpdp[2,0,2,0]= -(-3.00*G_1-24.00*G_3)
    U_dpdp[2,0,4,2]= -(-15.00*G_3)
    U_dpdp[3,0,0,2]= -(np.sqrt(3.00)*G_1-12.00*np.sqrt(3.00)*G_3)
    U_dpdp[3,0,1,2]= -(-3.00*G_1+6.00*G_3)
    U_dpdp[3,0,3,0]= -(-3.00*G_1-24.00*G_3)
    U_dpdp[3,0,4,1]= -(-15.00*G_3)
    U_dpdp[4,0,2,2]= -(-3.00*G_1+6.00*G_3)
    U_dpdp[4,0,3,1]= -(-3.00*G_1+6.00*G_3)
    U_dpdp[4,0,4,0]= -(-15.00*G_3)
    U_dpdp[0,1,0,1]= -(-G_1-18.00*G_3)
    U_dpdp[0,1,1,1]= -(-np.sqrt(3.00)*G_1-3.00*np.sqrt(3.00)*G_3)
    U_dpdp[0,1,2,0]= -(np.sqrt(3.00)*G_1+3.00*np.sqrt(3.00)*G_3)
    U_dpdp[0,1,4,2]= -(-2.00*np.sqrt(3.00)*G_1+9.00*np.sqrt(3.00)*G_3)
    U_dpdp[1,1,0,1]= -(-np.sqrt(3.00)*G_1-3.00*np.sqrt(3.00)*G_3)
    U_dpdp[1,1,1,1]= -(-3.00*G_1-24.00*G_3)
    U_dpdp[1,1,2,0]= -(-3.00*G_1+21.00*G_3)
    U_dpdp[1,1,4,2]= -(15.00*G_3)
    U_dpdp[2,1,0,0]= -(np.sqrt(3.00)*G_1+3.00*np.sqrt(3.00)*G_3)
    U_dpdp[2,1,1,0]= -(3.00*G_1-21.00*G_3)
    U_dpdp[2,1,2,1]= -(-3.00*G_1-24.00*G_3)
    U_dpdp[2,1,3,2]= -(-15.00*G_3)
    U_dpdp[3,1,2,2]= -(-3.00*G_1+6.00*G_3)
    U_dpdp[3,1,3,1]= -(-15.00*G_3)
    U_dpdp[3,1,4,0]= -(-3.00*G_1+6.00*G_3)
    U_dpdp[4,1,0,2]= -(np.sqrt(3.00)*G_1-12.00*np.sqrt(3.00)*G_3)
    U_dpdp[4,1,1,2]= -(3.00*G_1-6.00*G_3)
    U_dpdp[4,1,3,0]= -(-15.00*G_3)
    U_dpdp[4,1,4,1]= -(-3.00*G_1-24.00*G_3)
    U_dpdp[0,2,0,2]= -(-4.00*G_1-27.00*G_3)
    U_dpdp[0,2,3,0]= -(np.sqrt(3.00)*G_1-12.00*np.sqrt(3.00)*G_3)
    U_dpdp[0,2,4,1]= -(np.sqrt(3.00)*G_1-12.00*np.sqrt(3.00)*G_3)
    U_dpdp[1,2,1,2]= -(-15.00*G_3)
    U_dpdp[1,2,3,0]=-(-3.00*G_1+6.00*G_3)
    U_dpdp[1,2,4,1]=-(3.00*G_1-6.00*G_3)
    U_dpdp[2,2,2,2]= -(-15.00*G_3)
    U_dpdp[2,2,3,1]= -(-3.00*G_1+6.00*G_3)
    U_dpdp[2,2,4,0]= -(-3.00*G_1+6.00*G_3)
    U_dpdp[3,2,0,0]= -(-2.00*np.sqrt(3.00)*G_1+9.00*np.sqrt(3.00)*G_3)
    U_dpdp[3,2,1,0]= -(-15.00*G_3)
    U_dpdp[3,2,2,1]= -(-15.00*G_3)
    U_dpdp[3,2,3,2]= -(-3.00*G_1-24.00*G_3)
    U_dpdp[4,2,0,1]= -(-2.00*np.sqrt(3.00)*G_1+9.00*np.sqrt(3.00)*G_3)
    U_dpdp[4,2,1,1]= -(15.00*G_3)
    U_dpdp[4,2,2,0]= -(-15.00*G_3)
    U_dpdp[4,2,4,2]= -(-3.00*G_1-24.00*G_3)
    
    #for i in range(5):
    #    for j in range(5):
    #        for k in range(5):
    #            for l in range(5):
    #                U_multiplet[i,j,k,l] = U_rest[i,j,k,l] * 0.5
    #        for k in range(3):
    #            for l in range(3):
    #                U_multiplet[k,i,j,l] += U_pddp[k,i,j,l] 
    #                U_multiplet[i,k,j,l] += U_dpdp[i,k,j,l]
    ''' 
    do pp=1,5
    do mm=1,5
    do kk=1,5
    do jj=1,5
          U_multiplet(pp,mm,kk,jj) = U_rest(pp,mm,kk,jj)*0.5
    enddo
    enddo
    enddo
    enddo
    
    do pp=1,3
    do mm=1,5
    do kk=1,5
    do jj=1,3
          U_multiplet(pp+14,mm,kk,jj+14) = U_pddp(pp,mm,kk,jj)
          U_multiplet(mm,pp+14,kk,jj+14) = U_dpdp(mm,pp,kk,jj)
    enddo
    enddo
    enddo
    enddo
    '''

    #print 'U_ext', U_ext
    #print 'J_ext', J_ext
    #print 'U_pddp', U_pddp
    #print 'U_dpdp', U_dpdp

    return hop, U_onsite, U_ext, J_ext, U_rest, U_pddp, U_dpdp
