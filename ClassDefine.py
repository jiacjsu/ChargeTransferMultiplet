
# orders maps the order of the hr.dat into the orders of orbitals in this program
# example of use: 
#    for ii in range(N):
#        for jj in range(N):
#            hop[orders[ii] - 1, orders[jj] - 1] = value


class SpectraParas:
    def setSpectraParas(self, start_Eloss, end_Eloss, div_Eloss, \
                       start_Ein, end_Ein, div_Ein, \
                       epsilone_Ein, epsilone_Eloss, \
                       theta_in, phi_in, theta_out, phi_out):
        self.start_Eloss = start_Eloss
        self.end_Eloss = end_Eloss
        self.div_Eloss = div_Eloss
        self.start_Ein = start_Ein
        self.end_Ein = end_Ein
        self.div_Ein = div_Ein
        self.epsilone_Ein = epsilone_Ein
        self.epsilone_Eloss = epsilone_Eloss
        self.theta_in = theta_in
        self.phi_in = phi_in
        self.theta_out = theta_out
        self.phi_out = phi_out

        

class ModelParas:
   
    def setParas(self, F0, F2, F4, F0_pd, F2_pd, G1_pd, G3_pd, N_ligand, E_core, \
                 lambda_SO, seedname, orders, nup, ndn):
        self.rA = F0 - 49*F4/441
        self.rB = F2/49 - 5*F4/441
        self.rC = 35*F4/441
        self.F_0 = F0_pd
        self.F_2 = F2_pd/35
        self.G_1 = G1_pd/15
        self.G_3 = G3_pd/245
        self.N_ligand = N_ligand
        self.E_core = E_core
        self.lambda_SO = lambda_SO
        self.seedname = seedname
        self.orders = orders
        self.nup = nup
        self.ndn = ndn
    
    def printParas(self):
        print ''
        print '    rA =', self.rA, ', rB =', self.rB, ', rC =', self.rC, '(eV)'
        print '    F_0 =', self.F_0, ', F_2 =', self.F_2, ', G_1 =', \
            self.G_1, ', G_3 =', self.G_3
        print '    N_ligand =', self.N_ligand, ', E_core =', self.E_core, \
            ', lambda_SO =', self.lambda_SO, 
        print '    orders =', self.orders
        print '    nup =', self.nup, ', ndn =', self.ndn, ', seedname =', self.seedname
        print ''
        
    def get_rA(self):
        return self.rA
        
    def get_rB(self):
        return self.rB
    
    def get_rC(self):
        return self.rC   
        
    def get_F_0(self):
        return self.F_0
        
    def get_F_2(self):
        return self.F_2
        
    def get_G_1(self):
        return self.G_1
        
    def get_G_3(self):
        return self.G_3
        
    def get_N_ligand(self):
        return self.N_ligand
        
    def get_E_core(self):
        return self.E_core
                
    def get_lambda_SO(self):
        return self.lambda_SO
        
    def get_seedname(self):
        return self.seedname
        
    def get_orders(self):
        return self.orders
        
    def get_nup(self):
        return self.nup
        
    def get_ndn(self):
        return self.ndn
        
