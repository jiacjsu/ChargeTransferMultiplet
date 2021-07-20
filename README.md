# CTHFAM

Author: 
	Chunjing Jia

Date:
	June 6, 2018

What it does:

	This is the 64-bit state code for calculating the charge-transfer multiplet calculations.
	The on-site energys (includes ligand field splitting) and charge-transfer energy are 
	calculated from DFT+Wannier90. 
	Input files:
		seedname_hr.dat

How to compile:

	the C code
	gcc -shared -fPIC GenMatrixCore.c -o GenMatrixCore.so
	gcc -shared -fPIC GenHspCore.c -o GenHspCore.so

	Cython compile
	python setup.py build_ext --inplace

Running:

	python xasGUI.py

Examples:

LaNiO3

	9 ligands
	 total: 5 + 9 + 3 = 17 orbitals

	1  Ni dz2
	2  Ni dxz
	3  Ni dyz
	4  Ni dx2-y2
	5  Ni dxy

	25 valence electrons (39-14)

	3 + 12 = 15 ups, 3 + 13 = 16 spin down electrons




Working log
Feb 15 2018:
	Change based on src_Udholelanguage blablabla

	use ctypes for C function
	then call C function within Python


June 6 2018:
	Change unsigned int to unsigned long int, since when we have 9 ligands,
	we exceed the 32-bit int limit.

		When changing it to unsigned long int, please note that the intermediate results
		such as 1 << N (when N > 32) during the calculation, may run errors, because
		I assume that the intermediate states are saved only via regular int.
	(fixed)


		problem: ground state Hamiltonian w/o core-hole should not have spin-flip.
		but we see that, so there is a bug in the code...
	(fixed) ntmax = ntolmax
	
		
		problem: RIXS result does not converge.
		
		problem: When the Hilbert space is too small, lapack does not work# ChargeTransferMultiplet
