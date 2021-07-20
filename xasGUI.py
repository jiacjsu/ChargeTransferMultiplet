#!/usr/bin/python

# staticbox.py
# div_Ein format = float

import wx
from ClassDefine import ModelParas, SpectraParas
from Spectra import XAS, RIXS2
import numpy as np
import plotXAS_GUI as pxXAS
import plotRIXS_GUI as pxRIXS

class MyDialog(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(565, 660))



        # Section 1 Material Info

        wx.StaticBox(self, -1, 'Material Info', (10, 10), size=(270, 560))
        wx.StaticText(self, -1, 'Seedname', (25, 35))
        self.seedname = wx.TextCtrl(self, -1, '', (130, 32), size=(115, 25))
        wx.StaticText(self, -1, 'Transition Metal', (25, 70))
        self.metalname = wx.TextCtrl(self, -1, '', (160, 67), size=(85, 25))
        wx.StaticText(self, -1, 'dz2    dx2-y2   dxy      dxz       dyz', (25, 100))
        self.ind_dz2 = wx.TextCtrl(self, -1, '1', (25, 120), size=(30, 25))
        self.ind_dx2_y2 = wx.TextCtrl(self, -1, '2', (70, 120), size=(30, 25))
        self.ind_dxy = wx.TextCtrl(self, -1, '3', (115, 120), size=(30, 25))
        self.ind_dxz = wx.TextCtrl(self, -1, '4', (160, 120), size=(30, 25))
        self.ind_dyz = wx.TextCtrl(self, -1, '5', (205, 120), size=(30, 25))
        # may be changed into a xiala menu

        wx.StaticText(self, -1, 'Number of Ligand orbitals', (25, 160))
        self.N_ligand = wx.TextCtrl(self, -1, '', (190, 157), size=(55, 25))

        wx.StaticText(self, -1, 'Total Number of electrons', (25, 190))
        wx.StaticText(self, -1, 'Spin up', (25, 210))
        self.nup = wx.TextCtrl(self, -1, '', (75, 207), size=(45, 25))
        wx.StaticText(self, -1, 'Spin down', (130, 210))
        self.ndn = wx.TextCtrl(self, -1, '', (200, 207), size=(45, 25))

        #MaterialButton = wx.Button(self, -1, "Load Materials Info", (15, 225))       
        #MaterialButton.Bind(wx.EVT_BUTTON, self.loadMaterial)



        # Section 2 Material Parameters

        #wx.StaticBox(self, -1, 'Materials Parameters', (5, 265), size=(250, 335))
        self.loadAtomicButton = wx.Button(self, -1, "Load Atomic Values", (15, 252))       
        wx.StaticText(self, -1, 'with factor', (156, 255))
        self.factor = wx.TextCtrl(self, -1, '0.8', (226, 252), size=(35, 25))
        self.loadAtomicButton.Bind(wx.EVT_BUTTON, self.loadAtomicRacheParas)

        self.calculateWannierButton = wx.Button(self, -1, "Calculate from Wannier Orbitals", (15, 285))       


        wx.StaticText(self, -1, 'Rache Parameters for TM d orbitals', (25, 325))
        wx.StaticText(self, -1, 'F0_dd', (25, 345))
        self.F0 = wx.TextCtrl(self, -1, '', (40, 342), size=(55, 25))
        wx.StaticText(self, -1, 'F2_dd', (100, 345))
        self.F2 = wx.TextCtrl(self, -1, '', (115, 342), size=(55, 25))
        wx.StaticText(self, -1, 'F4_dd', (175, 345))
        self.F4 = wx.TextCtrl(self, -1, '', (190, 342), size=(55, 25))
   
        wx.StaticText(self, -1, 'Slater Integrals for TM dp interactions', (25, 380))
        wx.StaticText(self, -1, 'F0_pd', (25, 400))
        self.F0_pd = wx.TextCtrl(self, -1, '', (50, 397), size=(60, 25))
        wx.StaticText(self, -1, 'F2_pd', (125, 400))
        self.F2_pd = wx.TextCtrl(self, -1, '', (150, 397), size=(60, 25))
        wx.StaticText(self, -1, 'G1_pd', (25, 430))
        self.G1_pd = wx.TextCtrl(self, -1, '', (50, 427), size=(60, 25))
        wx.StaticText(self, -1, 'G3_pd', (125, 430))
        self.G3_pd = wx.TextCtrl(self, -1, '', (150, 427), size=(60, 25))

        
        wx.StaticText(self, -1, 'Core Energy Relative to Ef (eV)', (25, 470))
        self.E_core = wx.TextCtrl(self, -1, '', (220, 467), size=(40, 25))
        #wx.StaticText(self, -1, 'Core Hole Potential (eV)', (25, 505))
        #self.U_c = wx.TextCtrl(self, -1, '', (205, 502), size=(55, 25))
        wx.StaticText(self, -1, 'Spin-orbit Coupling lambda (eV)', (25, 540))
        self.lambda_SO = wx.TextCtrl(self, -1, '', (225, 537), size=(35, 25))

        self.MaterialButton = wx.Button(self, -1, "Load Material Parameters", (45, 575))       
        self.MaterialButton.Bind(wx.EVT_BUTTON, self.loadMaterial)
        


        # Section 3 Spectrum Info

        wx.StaticBox(self, -1, 'Spectrum Info', (290, 10), size=(265, 290))

        wx.StaticText(self, -1, 'Incoming X-ray Energies (eV)', (310, 35))
        wx.StaticText(self, -1, 'Start', (310, 55))
        self.start_Ein = wx.TextCtrl(self, -1, '', (345, 52), size=(40, 27))
        wx.StaticText(self, -1, 'End', (390, 55))
        self.end_Ein = wx.TextCtrl(self, -1, '', (420, 52), size=(40, 27))
        wx.StaticText(self, -1, 'Step', (465, 55))
        self.step_Ein = wx.TextCtrl(self, -1, '0.1', (500, 52), size=(40, 25))
        wx.StaticText(self, -1, 'Lorentzian Broadening (eV)', (310, 85))
        self.epsilone_Ein = wx.TextCtrl(self, -1, '0.2', (480, 82), size=(60, 25))

        wx.StaticText(self, -1, 'Energy Loss (eV)', (310, 120))
        wx.StaticText(self, -1, 'Start', (310, 140))
        self.start_Eloss = wx.TextCtrl(self, -1, '0.0', (345, 137), size=(40, 27))
        wx.StaticText(self, -1, 'End', (390, 140))
        self.end_Eloss = wx.TextCtrl(self, -1, '20.0', (420, 137), size=(40, 27))
        wx.StaticText(self, -1, 'Step', (465, 140))
        self.step_Eloss = wx.TextCtrl(self, -1, '0.05', (500, 137), size=(40, 25))
        wx.StaticText(self, -1, 'Lorentzian Broadening (eV)', (310, 170))
        self.epsilone_Eloss = wx.TextCtrl(self, -1, '0.05', (480, 167), size=(60, 25))

        wx.StaticText(self, -1, 'Incoming X-ray Polarizations', (310, 205))
        self.rb1 = wx.RadioButton(self, -1, 'x', (330, 225), style=wx.RB_GROUP)
        self.rb2 = wx.RadioButton(self, -1, 'y', (370, 225))
        self.rb3 = wx.RadioButton(self, -1, 'z', (410, 225))
        self.rb4 = wx.RadioButton(self, -1, 'all', (450, 225))

        wx.StaticText(self, -1, 'Outgoing X-ray Polarizations', (310, 255))
        self.rb5 = wx.RadioButton(self, -1, 'x', (330, 275), style=wx.RB_GROUP)
        self.rb6 = wx.RadioButton(self, -1, 'y', (370, 275))
        self.rb7 = wx.RadioButton(self, -1, 'z', (410, 275))
        self.rb8 = wx.RadioButton(self, -1, 'all', (450, 275))

        SpecButton = wx.Button(self, -1, "Load Spectrum Parameters", (325, 300))       
        SpecButton.Bind(wx.EVT_BUTTON, self.loadSpecParas)



        XASButton = wx.Button(self, -1, "Calculate L-edge XAS", (340, 360))       
        XASButton.Bind(wx.EVT_BUTTON, self.calculateXAS)
        RIXSButton = wx.Button(self, -1, "Calculate L-edge RIXS", (335, 395))       
        RIXSButton.Bind(wx.EVT_BUTTON, self.calculateRIXS)



        wx.Button(self, 1, 'Ok', (220, 605), (60, -1))

        self.Bind(wx.EVT_BUTTON, self.OnClose, id=1)

        self.ShowModal()
        self.Destroy()

    def OnClose(self, event):
        self.Close()

    def loadAtomicRacheParas(self, event):
        metalname = self.metalname.GetValue()
        if metalname == 'Fe':
            #self.rA.SetValue('1.8')
            #self.rB.SetValue('0.1258')
            #self.rC.SetValue('0.4652')
            #self.F_0.SetValue('0.0')
            #self.F_2.SetValue('0.155')
            #self.G_1.SetValue('0.26688')
            #self.G_3.SetValue('0.00928653')
            self.F0.SetValue('1.8')
            self.F2.SetValue('0.1258')
            self.F4.SetValue('0.4652')
            self.F0_pd.SetValue('0.0')
            self.F2_pd.SetValue('0.155')
            self.G1_pd.SetValue('0.26688')
            self.G3_pd.SetValue('0.00928653')
            self.E_core.SetValue('-705.0')
            self.lambda_SO.SetValue('8.5')
        if metalname == 'Ni':
            #self.rA.SetValue('2.0')
            #self.rB.SetValue('0.1517')
            #self.rC.SetValue('0.5294')
            #self.F_0.SetValue('0.0')
            #self.F_2.SetValue('0.19057')
            #self.G_1.SetValue('0.328')
            #self.G_3.SetValue('0.01142')
            self.F0.SetValue('0.5719')
            self.F2.SetValue('11.142')
            self.F4.SetValue('6.874')
            self.F0_pd.SetValue('0.448')
            self.F2_pd.SetValue('6.667')
            self.G1_pd.SetValue('4.922')
            self.G3_pd.SetValue('2.796')
            self.E_core.SetValue('-845.0')
            self.lambda_SO.SetValue('12.51')

            # Parameters from Haverkort's PRB

            # F2 = 11.14
            # F4 = 6.87
            # F2_2p3d = 6.67
            # G1_2p3d = 4.92
            # G3_2p3d = 2.80
            # lamda_2p = 11.51
            # lamda_3p = 1.4
            
            # U_3d,3d = 7.3
            # U_2p,3d = 8.5

   
            # Ni(2+) from the Journal of Physical Chemistry B, 2013, 117, 16512 by Kristjan Kunnus 
            # rC=0.42
            # rB=0.13
            # rA=2.0
            # Ni(2+) 3d-2p interaction from the about PRB
            # F_0=0.0d0; 
            # F_2=7.721/35.0*0.8; 
            # G_1=5.787/15.0*0.8; 
            # G_3=3.291/245.0*0.8;
            # lambda = 20 * 0.667


            # Translation relations
            # F_2=F^2/35 G_1=G^1/15; G_3=G^3/245;
            # F0 = 1/5 * (5rA + 7rC)
            # F2 = 1/7 * (7rB + rC) * 49
            # F4 = 21 * 21 * rC / 35
            # or rC = F4 * 35 / 21 / 21

    def loadMaterial(self, event):
        seedname =  self.seedname.GetValue()
	ind_dz2 = int(self.ind_dz2.GetValue())
	ind_dx2_y2 = int(self.ind_dx2_y2.GetValue())
	ind_dxy = int(self.ind_dxy.GetValue())
	ind_dxz = int(self.ind_dxz.GetValue())
	ind_dyz = int(self.ind_dyz.GetValue())
        N_ligand = int(self.N_ligand.GetValue())
        nup = int(self.nup.GetValue())
        ndn = int(self.ndn.GetValue())

        fac = float(self.factor.GetValue())
        F0 = float(self.F0.GetValue())
        F2 = float(self.F2.GetValue()) * fac
        F4 = float(self.F4.GetValue()) * fac
        F0_pd = float(self.F0_pd.GetValue())
        F2_pd = float(self.F2_pd.GetValue()) * fac
        G1_pd = float(self.G1_pd.GetValue()) * fac
        G3_pd = float(self.G3_pd.GetValue()) * fac
        E_core = -1*float(self.E_core.GetValue())
        #U_c = -1*float(self.U_c.GetValue())
        lambda_SO = float(self.lambda_SO.GetValue())
        
        orders0 = [ind_dz2, ind_dx2_y2, ind_dxy, ind_dxz, ind_dyz]
        orders = [0] * (5+N_ligand)
        for i in range(5):
            orders[orders0[i]-1] = i+1
        tag = 6
        for i in range(5+N_ligand):
            if orders[i] == 0:
                orders[i] = tag
                tag += 1

        self.para = ModelParas()
	#para.setParas(rA = 1.5, rB = 0.1365, rC = 0.5093, F_0 = 0.0, F_2 = 0.155, \
        #          G_1 = 0.26688, G_3 = 0.00928653, N_ligand = 3, E_core = 705.0, \
        #          lambda_SO = 8.5, U_c = 1.8, seedname = 'Ferricyanide', \
        #          orders = [5, 6, 1, 7, 4, 3, 2, 8], \
        #          nup = 6, ndn = 5)
        self.para.setParas(F0 = F0, F2 = F2, F4 = F4, F0_pd = F0_pd, F2_pd = F2_pd, \
                  G1_pd = G1_pd, G3_pd = G3_pd, N_ligand = N_ligand, E_core = E_core, \
                  lambda_SO = lambda_SO, seedname = seedname, \
                  orders = orders, \
                  nup = nup, ndn = ndn)
    

    def loadSpecParas(self, event):

        start_Ein = float(self.start_Ein.GetValue())
        end_Ein = float(self.end_Ein.GetValue())
        step_Ein = float(self.step_Ein.GetValue())
        epsilone_Ein = float(self.epsilone_Ein.GetValue())
        start_Eloss = float(self.start_Eloss.GetValue())
        end_Eloss = float(self.end_Eloss.GetValue())
        step_Eloss = float(self.step_Eloss.GetValue())
        epsilone_Eloss = float(self.epsilone_Eloss.GetValue())

        self.specPara = SpectraParas()
        #specPara.setSpectraParas(start_Eloss = 0.0, end_Eloss = 20.0, div_Eloss = 100, \
        #          start_Ein = 670.0, end_Ein = 730.0, div_Ein = 600, \
        #          epsilone_Ein = 0.2, epsilone_Eloss = 0.05, \
        #          theta_in = [np.pi/2.0, 0.0, 0.0], phi_in = [0.0, np.pi/2.0, 0.0], \
        #          theta_out = None, phi_out = None)
        self.specPara.setSpectraParas(start_Eloss = start_Eloss, end_Eloss = end_Eloss, \
                  div_Eloss = (end_Eloss - start_Eloss)/step_Eloss, \
                  start_Ein = start_Ein, end_Ein = end_Ein, div_Ein = (end_Ein - start_Ein)/step_Ein, \
                  epsilone_Ein = epsilone_Ein, epsilone_Eloss = epsilone_Eloss, \
                  theta_in = [np.pi/2.0, 0.0, 0.0], phi_in = [0.0, np.pi/2.0, 0.0], \
                  theta_out = [np.pi/2.0, 0.0, 0.0], phi_out = [0.0, np.pi/2.0, 0.0])
                  #theta_in = [0.0], phi_in = [0.0], \
                  #theta_out = [0.0], phi_out = [0.0])

    
    def calculateXAS(self, event):
        XAS(self.para, self.specPara)
        #app = wx.App(1)
        seedname = self.seedname.GetValue()
        app.frame = pxXAS.GraphFrame(seedname)
        app.frame.Show()
        app.MainLoop()

    def calculateRIXS(self, event):
        RIXS2(self.para, self.specPara)
        #app = wx.App(1)
        seedname = self.seedname.GetValue()
        app.frame = pxRIXS.GraphFrame(seedname)
        app.frame.Show()
        app.MainLoop()

if __name__ == '__main__':
    app = wx.App(0)
    MyDialog(None, -1, 'L-edge Spectrum calculator')
    app.MainLoop()
