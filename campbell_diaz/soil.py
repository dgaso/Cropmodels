# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 16:20:52 2019

@author: gaso001
"""
import sys
from pcse.base import SimulationObject, ParamTemplate, RatesTemplate, StatesTemplate, VariableKiosk
from pcse.decorators import prepare_rates, prepare_states
from pcse.traitlets import Float,Int, Instance, Enum, Unicode
import math
from mpmath import mp
from array import array
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import pandas as pd

from pcse.util import AfgenTrait


def convert_cm_to_m(v):
    return v/100.
def convert_KPa_to_hPa(p):
    return p*10.
def convert_hPa_to_KPa(k):
    return k/10.
def convert_ha_m2(m):
    return m*10000
def convert_j_Mj(j):
    return j/1000000
def convert_g_kg(g):
    return g/1000

    
class Water_balance(SimulationObject):
    '''Parameters**

    ============ ================================================= ==== ========
     Name        Description                                             Unit
    ============ ================================================= ==== ========
    FCP          Field capacity                                     mH2O m Soil
    PWPP         Permanent wilting point                            mH2O m Soi   
    ADWCP         Air dry water content                              mH2O m Soi   
    TCK          Thickness of soil layer                                    m
    RUNOFF1      Parameter 1 for runoff function   
    RUNOFF2      Parameter 2 for runoff function          
    RMIN         Root resistance                                         J kg-1 
    PSIPWP       Permanent wilting water potential                       J kg-1 
    PSIFC        Field capacity water potential                          J kg-1 
    S            Surface storage condition                                  m
    RDMAX        Maximun root depth                                         m
                                           m               
    ============ ================================================= ==== ========
    Rates**

    ============ ================================================= ==== ========
     Name        Description                                             Unit
    ============ ================================================= ==== ========
    FC           Parameter FCP in m of H2O per layer                       m
    PWP          Parameter PWPP in m of H2O per layer                      m
    INTERC       Precipitation intercepted by the canopy                    m
    RUNOFF       Rate of runoff                                             m
    INFIL        Water available for infiltration                           m
    RW           Rate of recharging water in each layer                     m
    NWC          Rate of water in the 1st layer                             m
    EVS          Rate of soil evaporation                                   m
    RDr          Rate of root growth                                        m
    B            Value for power equation of soil potential
    A            Value for power equation of soil potential
    SPSI         Soil potential per layer                            J kg-1 m-1
    FROOT        Root fraction per layer
    AVESPSI      Soil potential weighten by rooting fraction         J kg-1 m-1
    RBAR         Root resistance                                     J kg-1 m-1
    PSIX         Xilem Potential                                     J kg-1 m-1
    LOSS         Rate of water loss per layer                            mH2O m
    TL           Transpiration per layer                                    m
    AT           Actual Evapotranspiration per layer                        m
    T            Total Evapotranspiration                                   m
   
    ============ ================================================= ==== ========
    State variables**
 
    Name         Description                                             Unit
    ============ ================================================= ==== ========
    TINTERC      Total interception by the canopy                           m
    GPREC        Ground precipitation                                       m
    TRUNOFF      Total runoff                                               m
    PERC         Total deep percolation                                     m
    TWC          Water content                                          m-3 m-3
    TEVS         Total soil evaporation                                     m
    TRD          Root depth                                                 m
    TFR          Root fraction in each layer 
                                    
    ============ ================================================= ==== ========
    **External dependencies:**
    
    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
     FR        Fraction partitioned to roots.                     Y    -
     FS        Fraction partitioned to stems.                     Y    -
     FL        Fraction partitioned to leaves.                    Y    -
     FO        Fraction partitioned to storage orgains            Y    -
     DVS       Development stage
    =======  =================================== =================  ============    
    '''
    pheno = Instance(SimulationObject)

    class Parameters(ParamTemplate):                        
        #Soil water balance
        FCP    = Float(-99.)
        PWPP   = Float(-99.)
        ADWCP  = Float(-99.)
        TCK    = Float(-99.)
        RUNOFF1 = Float(-99.)
        RUNOFF2= Float(-99.)
        RMIN  = Float(-99.)
        PSIPWP= Float(-99.)
        PSIFC = Float(-99.)
        S     = Float(-99.)
        
        RDMAX = Float(-99.)
        
    class RateVariables(RatesTemplate):
        PE        = Float(-99.)
        PT        = Float(-99.)
        EPT        = Float(-99.)
        PREC      = Float(-99.)
        nl        = Int
        FC        = Float(-99.) 
        PWP       = Float(-99.) 
        ADWC      = Float(-99.)
        
        B         = Float(-99.)
        A         = Float(-99.)        
        INTERC    = Float(-99.)        
        GPREC     = Float(-99.)        
        RUNOFF    = Float(-99.)
        INFIL     = Float(-99.)        
        RW        = Float(-99.)        
        values_RW = [] 
        RWATER    = Instance(np.ndarray)         
        NWC       = Float(-99.)
        EVS       = Float(-99.)                       
        SPSI      = Float(-99.)
        values_SPSI = []
        SP        = Instance(np.ndarray) 
        #root
        FcR       = Instance(np.ndarray) 
        values_FR = []
        FROOT     = Float(-99.)
        AVEPSI    = Float(-99.)
        RBAR      = Float(-99.)
        PSIX      = Float(-99.) 
        z         = Float(-99.) 

        LOSS      = Instance(np.ndarray) 
        W         = Instance(np.ndarray)                                  
        arr_bool  = Instance(np.ndarray) 
        arr_TRUE  = Instance(np.ndarray) 
        TL        = Instance(np.ndarray) 
        net_RW    = Instance(np.ndarray) 
        T         = Float(-99.)
        arr_subt  = Instance(np.ndarray) 
        arr_mult  = Instance(np.ndarray) 
        AT         = Instance(np.ndarray) 
        check = Float(-99.)

    class StateVariables(StatesTemplate):  
        PTa     = Float(-99.)
        #Soil water balance
        TPE     = Float(-99.)
        TPT     = Float(-99.)
        TPREC   = Float(-99.)
        TINTERC = Float(-99.)
        TRUNOFF = Float(-99.)         
        PERC    = Float(-99.)                
        TE      = Float(-99.)
        Ta      = Float(-99.)    
        WC      = Instance(np.ndarray) 
        WCv     = Instance(np.ndarray) 
        #To chech WB closed
        TT        = Float(-99.)
        WB_close  = Float(-99.)
        Diff_WC   = Float(-99.)
        TWC       = Float(-99.)
        W_Stress  = Float(-99.)
                
    def initialize(self, day, kiosk, parametervalues):        
        self.params = self.Parameters(parametervalues)
        self.rates = self.RateVariables(kiosk,publish = None)        
        self.kiosk = kiosk
        layers = math.floor(self.params.RDMAX/self.params.TCK)
        self.states = self.StateVariables(kiosk, publish=["Ta", "W_Stress","PTa"],  PTa=0,
                                          TPE=0.0, TPT=0.0,
                                          WC=np.full(layers,self.params.FCP*self.params.TCK), 
                                          WCv=np.full(layers,self.params.FCP*self.params.TCK),
                                          TT=0.0, Diff_WC=0., W_Stress=0, TINTERC=0, TRUNOFF=0,
                                          PERC=0, TE=0, Ta=0,  WB_close =0, TPREC=0.0, TWC=0.)

    @prepare_rates
    def calc_rates(self, day, drv):                
        p = self.params
        r = self.rates
        s = self.states
        k = self.kiosk     
        
        if "FI" not in self.kiosk: 
             k.FI = 0.0
        else: 
            k.FI=self.kiosk["FI"]
         
        if "TRD" not in self.kiosk: 
            k.TRD = 0.0
        else:
            k.TRD = self.kiosk["TRD"]

        r.PE = convert_cm_to_m(drv.ET0 * (1 - k.FI))
        r.PT = convert_cm_to_m(drv.ET0 * k.FI)
        r.ETP=convert_cm_to_m(drv.ET0 )*100*10
        

        # Convert fc and pwp to m of H2O
        r.nl = math.floor(p.RDMAX/p.TCK) 
        r.FC = p.FCP * p.TCK
        r.PWP = p.PWPP * p.TCK
        r.ADWC = p.ADWCP * p.TCK
        
        # Parameters for eq os SPSI
        r.B = math.log(p.PSIPWP/p.PSIFC)/math.log(p.FCP/p.PWPP)
        r.A = math.exp((math.log(-p.PSIFC)) + r.B * math.log(max(p.FCP , sys.float_info.min)) )

        # Interception
        r.PREC = drv.RAIN /100
       
        if drv.RAIN != 0.:
            r.INTERC = min(drv.RAIN/100, 0.001 * k.FI) 
        r.GPREC = ((drv.RAIN/100) - r.INTERC)
        
        # Runoff calculation
        if r.GPREC <= 0.2 * p.S:
            r.RUNOFF = 0
        else:
            r.RUNOFF = ((r.GPREC - p.RUNOFF1 * p.S)**2)/(r.GPREC + p.RUNOFF2 * p.S)                
        r.INFIL = ((drv.RAIN/100) - r.INTERC - r.RUNOFF)

        # Recharge water per layer
        r.values_RW = []
        r.values_FR = []
        r.values_SPSI = []
        for j in range (s.WC.shape[0]):
            if r.INFIL>0:
                if r.INFIL <= (r.FC - s.WC[j]):
                    r.RW = r.INFIL
                    r.INFIL = 0.
                else:                 
                    r.RW = r.FC - s.WC[j] 
                    r.INFIL =  r.INFIL - (r.FC - s.WC[j])    
            else: 
                r.RW = 0                
                            
            r.values_RW.append(r.RW)         
            r.RWATER = np.array(r.values_RW)  
            
            # Evaporation calculation
            if j < 1:
                # if s.WC[0]>r.PWP:
                #      r.EVS=r.PE
                if s.WC[0] < r.PWP:
                    r.PE = r.PE * (((s.WC[0]/p.TCK)-p.ADWCP)/(p.PWPP - p.ADWCP))**2.
                r.NWC = s.WC[0] - r.PE  
                r.EVS = r.PE               
                if r.NWC < r.ADWC: 
                    r.NWC = r.ADWC                    
                    r.EVS = r.EVS - r.NWC
                                
            # Fraction root per layer
            if j >= 1:
                r.z+=p.TCK
                if r.z <= k.TRD: 
                   r.FROOT = p.TCK*(2.*(k.TRD - r.z) + p.TCK)/(k.TRD*k.TRD)
                elif r.z > k.TRD and (r.z - p.TCK) < k.TRD:
                   r.FROOT = ((k.TRD - r.z + p.TCK)/k.TRD)**2.
                else:
                    r.FROOT = 0.

                r.SPSI = -r.A*mp.exp(-r.B*math.log(max((s.WC[j]/p.TCK), sys.float_info.min)))
                r.AVEPSI += (r.FROOT*r.SPSI)                                      
                r.values_FR.append(r.FROOT)  
                r.FcR = np.array(r.values_FR)
                r.values_SPSI.append(r.SPSI)
                r.SP = np.array(r.values_SPSI)

        r.RBAR = p.RMIN/max(k.FI, 1e-70)
        r.PSIX = r.AVEPSI - (r.RBAR*r.PT) 
    
        if r.PSIX < p.PSIPWP: r.PSIX=p.PSIPWP
       
        # Transpiration calculation
        r.arr_subt = np.subtract(r.SP, r.PSIX)
        r.arr_mult = np.multiply(r.FcR, r.arr_subt)
        r.LOSS = np.true_divide(r.arr_mult, (r.RBAR*p.TCK))
     
        # check loss under water deficit
        r.W = np.true_divide(s.WC[1::], p.TCK)
        r.arr_bool = (np.subtract(r.W, r.LOSS)) < p.PWPP
        r.arr_TRUE = np.subtract(r.W, p.PWPP)
        
        # Actual evapotranspiration        
        r.TL= np.multiply(np.where(r.arr_bool, r.arr_TRUE, r.LOSS), p.TCK)
        r.AT = np.append(r.EVS, r.TL)
        r.T = np.sum(r.TL)
        r.net_RW = np.subtract(r.RWATER, r.AT)
                

    @prepare_states
    def integrate(self, day,delt):
        p = self.params
        r = self.rates
        s = self.states
        k = self.kiosk

        s.PTa = r.PT                                                      
        # Soil water balance
        s.TPE += r.PE
        s.TPT = r.PT                                        
        s.TPREC=r.PREC
        s.TINTERC +=  r.GPREC
        s.TRUNOFF +=  r.RUNOFF
        s.PERC +=r.INFIL        
        s.TE += r.EVS
        s.Ta = r.T       
        s.WC = np.add(s.WC, r.net_RW)
        s.WCv = (s.WC / p.TCK)
        
          
        # check WB closed
        s.TT+=r.T
        s.Diff_WC = np.sum(np.subtract(np.full(r.nl, r.FC), s.WC) )
        s.WB_close = s.TINTERC - s.TRUNOFF - s.PERC - s.TT - s.TE + s.Diff_WC
        s.TWC = (np.sum(np.subtract(s.WC, r.PWP)))*1000
        
        s.W_Stress = r.T / max(r.PT, 1e-70)
        
