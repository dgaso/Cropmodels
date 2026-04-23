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
    RUNOFF1      Parameter 1 for runoff function                        -
    RUNOFF2      Parameter 2 for runoff function                        -
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
    INTERC       Precipitation intercepted by the canopy                    m
    GPREC        Ground precipitation                                       m
    RUNOFF       Rate of runoff                                             m
    INFIL        Water available for infiltration                           m
    EVS          Rate of soil evaporation                                   m
    TL           Transpiration per layer                                    m
    AT           Actual Evapotranspiration per layer                        m
    T            Total Evapotranspiration                                   m
    net_RW	     Net water recharge                                         m
    ============ ================================================= ==== ========
    Intermediate variables within Rates function**

    ============ ================================================= ==== ========
     Name        Description                                             Unit
    ============ ================================================= ==== ========
    RW           Rate of recharging water in each layer                     m
    NWC          Rate of water in the 1st layer                             m
    RDr          Rate of root growth                                        m
    B            Value for power equation of soil potential
    A            Value for power equation of soil potential
    SPSI         Soil potential per layer                            J kg-1 m-1
    FROOT        Root fraction per layer                                   -
    AVESPSI      Soil potential weighten by rooting fraction         J kg-1 m-1
    RBAR         Root resistance                                     J kg-1 m-1
    PSIX         Xilem Potential                                     J kg-1 m-1
    LOSS         Rate of water loss per layer                            mH2O m
       
    ============ ================================================= ==== ========   

 ============ ================================================= ==== ========
    State variables**
 
    Name         Description                                             Unit
    ============ ================================================= ==== ========
    TINTERC      Total interception by the canopy                           m
    TRUNOFF      Total runoff                                               m
    PERC         Total deep percolation                                     m   
    TE           Total soil evaporation                                     m
    WC           Water content per layer                                    m
    WCv          Volumetric water content per layer                     mH2O m Soil
    TT           Total transpiration                                        m
    WB_close     Check water balance                                        m
    TWC          Water content                                          m-3 m-3
                                    
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
        INTERC    = Float(-99.)        
        GPREC     = Float(-99.)        
        RUNOFF    = Float(-99.)
        INFIL     = Float(-99.)        
        EVS       = Float(-99.)                           
        T         = Float(-99.)
        TL        = Instance(np.ndarray) 
        AT        = Instance(np.ndarray) 
        net_RW    = Instance(np.ndarray)
        
    class StateVariables(StatesTemplate):  
        PTa     = Float(-99.)
        Ta      = Float(-99.)   
        Tp      = Float(-99.)
        #Soil water balance
        WC      = Instance(np.ndarray) 
        WCv     = Instance(np.ndarray)                    
        TINTERC = Float(-99.)
        TRUNOFF = Float(-99.)         
        PERC    = Float(-99.)                      
        TT        = Float(-99.)
        TE        = Float(-99.)
        WB_close  = Float(-99.)
        TWC       = Float(-99.)
        W_Stress  = Float(-99.)
        
        # Cumulaive variables to extract features      
        CRainv = Float(-99.)
        CRainr = Float(-99.)
        TWCR1 = Float(-99.)
        TWCR5 = Float(-99.)
                
    def initialize(self, day, kiosk, parametervalues):        
        self.params = self.Parameters(parametervalues)         
        self.rates = self.RateVariables(kiosk,publish = ["T","PT"]) 
        self.kiosk = kiosk        
        layers = math.floor(self.params.RDMAX/self.params.TCK)
        self.states = self.StateVariables(kiosk, publish=["Ta", "W_Stress","PTa"], Ta=0, PTa=0, Tp=0,                                        
                                          WC=np.full(layers,self.params.FCP*self.params.TCK), 
                                          WCv=np.full(layers,self.params.FCP*self.params.TCK),
                                          TT=0.0, Diff_WC=0., W_Stress=1, TINTERC=0, TRUNOFF=0,
                                          PERC=0, TE=0, WB_close =0,  
                                          TWC=layers * ((self.params.FCP - self.params.PWPP) * self.params.TCK) * 1000.0,
                                          CRainv =0., CRainr = 0.,TWCR1 = None, TWCR5 = None)


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
                   
        # Convert fc and pwp to m of H2O
        
        FC = p.FCP * p.TCK
        PWP = p.PWPP * p.TCK
        ADWC = p.ADWCP * p.TCK
               
        # Interception
        if drv.RAIN != 0.:
            r.INTERC = min(drv.RAIN/100, 0.001 * k.FI) 
        r.GPREC = ((drv.RAIN/100) - r.INTERC)
        
        # Runoff calculation
        if r.GPREC <= 0.2 * p.S:
            r.RUNOFF = 0
        else:
            r.RUNOFF = ((r.GPREC - p.RUNOFF1 * p.S)**2)/(r.GPREC + p.RUNOFF2 * p.S)                
        r.INFIL = ((drv.RAIN/100) - r.INTERC - r.RUNOFF)

        # Parameters for eq os SPSI
        B = math.log(p.PSIPWP/p.PSIFC)/math.log(p.FCP/p.PWPP)
        A = math.exp((math.log(-p.PSIFC)) + B * math.log(max(p.FCP , sys.float_info.min)) )
        AVEPSI = 0  
        z=0     
        
        # Recharge water per layer
        values_RW = []
        values_FR = []
        values_SPSI = []              
        for j in range (s.WC.shape[0]):
            if r.INFIL>0:
                if r.INFIL <= (FC - s.WC[j]):
                    RW = r.INFIL
                    r.INFIL = 0.
                else:                 
                    RW = FC - s.WC[j] 
                    r.INFIL =  r.INFIL - (FC - s.WC[j])    
            else: 
                RW = 0                
                            
            values_RW.append(RW)         
            RWATER = np.array(values_RW)  
            
            # Evaporation calculation
            if j < 1:
                if s.WC[0]>PWP:
                    r.EVS = r.PE                  
                else:                    
                    evapo = r.PE * (((s.WC[0]/p.TCK)-p.ADWCP)/(p.PWPP - p.ADWCP))**2. 
                    r.EVS = min(evapo, 1e6) 
                                           
                NWC = s.WC[0] - r.EVS
                             
                if NWC < ADWC: 
                    NWC = ADWC                    
                    s.WC[0]= NWC
                                           
            # Fraction root per layer           
            if j >= 1:               
                z+=p.TCK
                if z <= k.TRD:                     
                    FROOT = p.TCK*(2.*(k.TRD - z) + p.TCK)/(k.TRD*k.TRD)
                elif z > k.TRD and (z - p.TCK) < k.TRD:
                    FROOT = ((k.TRD - z + p.TCK)/k.TRD)**2.
                else:
                    FROOT = 0.

                # Accumualte
                SPSI = -A*mp.exp(-B*math.log(max((s.WC[j]/p.TCK), sys.float_info.min)))    
                AVEPSI += (FROOT*SPSI)  
                
                values_FR.append(FROOT)  
                FcR = np.array(values_FR)
                values_SPSI.append(SPSI)
                SP = np.array(values_SPSI)

        
        if k.FI <= 0.0:           
            RBAR = np.inf
            PSIX = p.PSIPWP
            LOSS = np.zeros_like(FcR, dtype=float)
            r.TL = np.zeros_like(FcR, dtype=float)
        else:   
            #RBAR = p.RMIN/max(k.FI, 1e-70)
            RBAR = p.RMIN / k.FI
            PSIX = AVEPSI - (RBAR*r.PT) 
        
            if PSIX < p.PSIPWP: PSIX=p.PSIPWP
           
            # Transpiration (LOSS) calculation
            arr_subt = np.subtract(SP, PSIX)
            arr_mult = np.multiply(FcR, arr_subt)
            LOSS = np.true_divide(arr_mult, (RBAR*p.TCK))
         
            # check loss under water deficit
            W = np.true_divide(s.WC[1::], p.TCK) #convert water content in m to volumentric
            arr_bool = (np.subtract(W, LOSS)) < p.PWPP
            arr_TRUE = np.subtract(W, p.PWPP)
            
            # Actual evapotranspiration        
            r.TL= np.multiply(np.where(arr_bool, arr_TRUE, LOSS), p.TCK)
            r.AT = np.append(r.EVS, r.TL)
            r.T = np.sum(r.TL)
            r.net_RW = np.subtract(RWATER, r.AT)
            
           

    @prepare_states
    def integrate(self, day,delt):
        p = self.params
        r = self.rates
        s = self.states
        k = self.kiosk

        nl = math.floor(p.RDMAX/p.TCK)
        FC = p.FCP * p.TCK
        PWP = p.PWPP * p.TCK   
        
        s.PTa = r.PT    
        s.Ta = r.T
        s.Tp = r.PT
        
        # Soil water balance  
        s.TINTERC +=  r.GPREC
        s.TRUNOFF +=  r.RUNOFF
        s.PERC +=r.INFIL        
        s.TE += r.EVS      
        s.WC = np.add(s.WC, r.net_RW)
        s.WCv = (s.WC / p.TCK)    
          
        # check WB closed
        s.TT+=r.T
        Diff_WC = np.sum(np.subtract(np.full(nl, FC), s.WC) )
        s.WB_close = s.TINTERC - s.TRUNOFF - s.PERC - s.TT - s.TE + Diff_WC
        s.TWC = (np.sum(np.subtract(s.WC, PWP)))*1000
        
        if r.PT <= 0.0 or k.TRD <= 0.0:
            s.W_Stress = 1.0
        else:
            s.W_Stress = r.T / r.PT
            
        #s.W_Stress = r.T / max(r.PT, 1e-70)
        
        if k.DVS > 0 and k.DVS <= 1.:
            s.CRainv += r.GPREC * 1000
               
        if k.DVS > 1 and k.DVS <= 2.:
            s.CRainr += r.GPREC * 1000
       
        if s.TWCR1 is None and k.DVS >= 1.:
            s.TWCR1 = max(0,s.TWC)
            
        if s.TWCR5 is None and k.DVS >= 1.5:
            s.TWCR5 = max(0,s.TWC)
            
       
            
            
            
            
            
        
