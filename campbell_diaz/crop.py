# -*- coding: utf-8 -*-
# Copyright (c) 2020 Wageningen-UR
# Deborah Gaso Melgar, June 2020

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

from .partitioning import DVS_Partitioning as Partitioning
from .wofost_soybean_phenology import SoybeanPhenology as Phenology
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


class Campbell(SimulationObject):
    """Parameters**

    ============ ================================================= ==== ========
     Name        Description                                             Unit
    ============ ================================================= ==== ========
    DVV1         Paramter 1 for vapor deficit equation
    DVV2         Paramter 2 for vapor deficit equation
    DVV3         Paramter 3 for vapor deficit equation
    NTR          Nitrogen content in grain                                g.g-1
    LNTR         Nitrogen content in leaves                               g.g-1
    FNTR         Fraction of N translocated from leaves to seeds
    HD           Factor to standardize humidity content at 13% 
    K            Light extinction coefficient (Kukal and Irmak, 2020)
    Ppar         Proportion of PAR in the total radiation
    WUE          Dry matter water ratio                                     Pa
    DSLA         Specific leaf area for dead leaves                    m-2 kg-1
    RUE          Radiation use efficiency (Kukal and Irmak, 2020)       g Mj m2
    GCC          Conversion coefficient (CHO to soyeban seed)
    LAIC         Critic leaf area index                                 m-2 m-2
    FTRANSL      Fraction of TDM to be translocated
    RDRSHM       Maximum relative death rate of leaves due to shading (LINTUL2)  d-1
    ============ ================================================= ==== ========
    Rates**

    ============ ================================================= ==== ========
     Name        Description                                             Unit
    ============ ================================================= ==== ========
    FI           Fractional interception
    PE           Potential evaporation                                      m
    PT           Potential transpiration                                    m
    VDD          Correction by vapor deficit water                         KPa
    PARi         Intercepted PAR                                         Mj m-2
    DM           Rate of growth                                          kg m-2
    ROOT         Growth rate root                                        kg m-2
    STEMS        Growth rate stems                                       kg m-2
    LEAF         Growth rate leaf                                        kg m-2
    SEED         Growth rate storage organs                              kg m-2
    TN           Translocated nitrogen from leaves to grains             kg m-2
    DLEAF        Senescence rate of leaf                                 kg m-2
    RDRSH        Relative death rate of leaves due to shading (LINTUL2)   d-1
    RDRT         Table of RDR as a function of temperature
    ============ ================================================= ==== ========
    State variables**

    Name         Description                                             Unit
    ============ ================================================= ==== ========
    TPE          Total potential evaporation                                m
    TPT          Total potential transpiration                              m
    TDM          Total above-ground biomass                              kg m-2
    TROOT        Total weight of roots                                   kg m-2
    TSTEMS       Total weight of stems                                   kg m-2
    TLEAF        Total weight of leaves                                 m-2 m-2
    TSEED        Total weight of storage organs                          kg m-2
    LAI          Leaf area index                                        m-2 m-2
    TDLEAF       Total of dead leaves                                   m-2 m-2
    GLEAF        Total of green leaves                                  m-2 m-2
    SLA          Specific leaf area                                    m-2 kg-1

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
    """
    # sub-model components for crop simulation
    pheno = Instance(SimulationObject)
    part  = Instance(SimulationObject)
    soil  = Instance(SimulationObject)

    class Parameters(ParamTemplate):                        
        RDRT   = AfgenTrait()
        RDRSHM = Float(-99.)
        LAIC   = Float(-99.)
        
        initLAI = Float(-99.)

        K      = Float(-99.)
        Ppar   = Float(-99.)
        WUE    = Float(-99.)
       
        DSLA   = Float(-99.) 
        NTR    = Float(-99.)
        LNTR   = Float(-99.)
        FNTR  = Float(-99.)
        HD     = Float(-99.)

        GCC = Float(-99.)
        FTRANSL = Float(-99.)
        
        SLATB   = AfgenTrait()
        RUE    = Float(-99.)   
        
        RDMAX = Float(-99.)
        
                              
    class RateVariables(RatesTemplate):

        Es_mn   = Float(-99.)
        Es_mx   = Float(-99.)
        Es_avg  = Float(-99.)
        VPD     = Float(-99.)
        
        RDRDV = Float(-99.)
        RDRSH = Float(-99.)
        RDR   = Float(-99.)
        DLAI  = Float(-99.)
        
        DM_W    = Float(-99.)
        DM_R    = Float(-99.)
        DM      = Float(-99.)  
        PDM      = Float(-99.)             
       
        
        FI      = Float(-99.)
        ROOT    = Float(-99.) 
        STEMS   = Float(-99.) 
        LEAF    = Float(-99.) 
        WLEAF    = Float(-99.) 
        SEED    = Float(-99.)
        PSEED   = Float(-99.)
        TN      = Float(-99.) 
        WDLEAF  = Float(-99.)
        DLEAF   = Float(-99.)
        GLEAF   = Float(-99.)
        TRANSL  = Float(-99.)
        
        #root
        RD        = Float(-99.)
        WD        = Float(-99.)
       

    class StateVariables(StatesTemplate):
        TDM     = Float(-99.)
        TSTEM   = Float(-99.) 
        TLEAF   = Float(-99.)        
        TSEED   = Float(-99.) 
        YIELD   = Float(-99.) 
        LAI     = Float(-99.) 
        TDLEAF  = Float(-99.)        
        SLA     = Float(-99.)

        TDMTRANSL = Float(-99.)
        POOLTRSL = Float(-99.)
        #root
        TRD     = Float(-99.)
        da      = Int
        #cumulative variable to extract features
        LAIR1 = Float(-99.) 
        LAIR5 = Float(-99.) 
        TDMFlowering = Float(-99.)
        TDMR1 = Float(-99.)
        TDMR5 = Float(-99.)
        CWDv  = Float(-99.)
        CWDr  = Float(-99.)
        CVPDv  = Float(-99.)
        CVPDr  = Float(-99.)
        CTv  = Float(-99.)
        CTr  = Float(-99.)
        RADv = Float(-99.)
        RADr = Float(-99.)
       

        
    def initialize(self, day, kiosk, parametervalues):
        self.params = self.Parameters(parametervalues)
        self.rates = self.RateVariables(kiosk, publish=["FI"])        
        self.kiosk = kiosk                                                                                
        # Initialize components of the crop        
        self.pheno = Phenology(day, kiosk, parametervalues)
        self.part  = Partitioning(day, kiosk, parametervalues)
        
        DVS = self.kiosk["DVS"]
        SLA = self.params.SLATB(DVS)
        
        # =============================================================================
        #         # Initial total (living+dead) above-ground biomass of the crop
        #         FR = self.kiosk["FR"]
        #         FS = self.kiosk["FS"]
        #         FL = self.kiosk["FL"]
        #         FO = self.kiosk["FO"]
        # =============================================================================
        self.states = self.StateVariables(kiosk, publish=["TRD"],
                                       
                                          TDM=0.00,  GLEAF=0.0, TSTEM=0.0,
                                          TLEAF=0.0,TSEED=0.0,YIELD=0.0, 
                                          LAI=self.params.initLAI, TDLEAF =0, SLA=SLA, TRD=0.0,da=0,
                                          TDMFlowering=None, TDMR1=None, LAIR1=None,TDMR5=None, LAIR5=None,
                                           CVPDv = 0., CVPDr = 0., CTv= 0., CTr = 0., RADv = 0., RADr =0., 
                                          TDMTRANSL=0,POOLTRSL=0, CWDv=0.,CWDr=0.)

    @prepare_rates
    def calc_rates(self, day, drv):                
        p = self.params
        r = self.rates
        s = self.states
        k = self.kiosk
        
        self.pheno.calc_rates(day, drv)
        crop_stage = self.pheno.get_variable("STAGE")

        # if before emergence there is no need to continue
        # because only the phenology is running.
        if crop_stage == "emerging":
            return

        self.part.calc_rates(day, drv)
              
        r.Es_mn=0.6108 * math.exp( max(((17.27 * drv.TMIN)/(drv.TMIN + 237.3)),0.001) )
        r.Es_mx=0.6108 * math.exp( max(((17.27 * drv.TMAX)/(drv.TMAX + 237.3)),0.001) )
        r.Es_avg=(r.Es_mn + r.Es_mx )/2
        r.VPD=r.Es_avg - convert_hPa_to_KPa(drv.VAP)
            
        r.FI = 1. - mp.exp(-p.K * s.LAI)
        r.PARi = convert_j_Mj(drv.IRRAD) * p.Ppar * r.FI
                   
        if k.DVS < 2:
            if "Ta" not in self.kiosk:
                k.Ta = 0.001
            else:
                k.Ta = self.kiosk["Ta"]   
                
            if "W_Stress" not in self.kiosk:
                k.W_Stress = 0.0
            else:
                k.W_Stress = self.kiosk["W_Stress"]   
                
            if "PTa" not in self.kiosk:
                k.PTa = 0.0
            else:
                k.PTa = self.kiosk["PTa"]     
                             
            r.DM_W = k.Ta * (p.WUE/r.VPD)
            r.DM_R = convert_g_kg(r.PARi * p.RUE )
           
            r.DM = min(r.DM_W, r.DM_R) 
            
            r.PDM = k.PTa * (p.WUE/r.VPD)                  
            r.STEMS = r.DM * k.FS             
            r.WLEAF = r.DM * k.FL
            r.LEAF = r.DM * k.FL * convert_ha_m2(s.SLA)            
            r.PSEED = r.PDM * k.FO  
       
            # Biomass reallocated from vegetative to seed
            if s.TDMTRANSL>0:
                r.TRANSL = r.PSEED - (r.DM * k.FO)
                r.SEED = (r.DM * k.FO) + r.TRANSL
            else: r.SEED = r.DM * k.FO
            
            # Senescence from N translocation
            r.TN = ((r.SEED*p.NTR)/p.FNTR)  
            
             #senescence rate from LINTUL
            if k.DVS>1.5:              
                r.RDRDV = p.RDRT(drv.TEMP)
                r.RDRSH = p.RDRSHM *((s.LAI - p.LAIC)/ p.LAIC)
                if r.RDRSH >0 or r.RDRSH <0.03:
                    r.RDRSH = r.RDRSH 
                if r.RDRSH < 0:
                    r.RDRSH =0
                if r.RDRSH >0.03:
                    r.RDRSH =0.03          
            else: r.RDRDV = 0
                            
            r.RDR = max(r.RDRDV, r.RDRSH)   
            r.DLAI = s.LAI * r.RDR                                     
            r.WDLEAF = r.TN/p.LNTR                   
            r.DLEAF = r.WDLEAF * p.DSLA                                                      
            r.GLEAF = r.LEAF - max(r.DLAI, r.DLEAF)

        #Rooting growth
        r.RD = p.RDMAX * (1./(1+44.2*math.exp(-15*(s.da)/(140))))  
        r.WD = k.PTa - k.Ta
    

    @prepare_states
    def integrate(self,day,delt):
        p = self.params
        r = self.rates
        s = self.states
        k = self.kiosk

        # crop stage before integration
        crop_stage = self.pheno.get_variable("STAGE")        
        self.pheno.integrate(day, delt=1.0)
        # if before emergence there is no need to continue
        # because only the phenology is running.
        # Just run a touch() to to ensure that all state variables are available
        # in the kiosk
        if crop_stage == "emerging":
            self.touch()
            return

        self.part.integrate(day, delt=1.0)

        DVS = self.kiosk["DVS"]
        s.SLA = p.SLATB(DVS)  
   
        s.TDM += r.DM
        
        if s.TDMFlowering is None and k.DVS >= 1.:
            s.TDMFlowering = s.TDM
            s.TDMTRANSL = s.TDMFlowering * p.FTRANSL            
            s.TDMR1 = s.TDM * 10000
               
        s.TDMTRANSL -= r.TRANSL
        s.TSEED += r.SEED         
        s.YIELD = s.TSEED * p.GCC * p.HD * 10000
        s.TSTEM += r.STEMS
        s.TLEAF += r.DM * k.FL 
        s.LAI +=  r.GLEAF
        #print(day,k.DVS)
        s.da+=1  
        s.TRD = r.RD
        
        ##Cumulative variable to extract features
                                      
        if s.LAIR1 is None and k.DVS >= 1.:
            s.LAIR1 = s.LAI            
        if s.LAIR5 is None and k.DVS >= 1.5:
            s.LAIR5= s.LAI

        if s.TDMR5 is None and k.DVS >= 1.5:
            s.TDMR5 = s.TDM * 10000
        
        if k.DVS > 0 and k.DVS <= 1.:
            s.CVPDv += r.VPD 
            s.CWDv += r.WD * 1000
            s.CTv += k.Ta * 1000
            s.RADv += r.PARi
            
        if k.DVS > 1 and k.DVS <= 2.:
            s.CVPDr += r.VPD
            s.CWDr += r.WD * 1000
            s.CTr += k.Ta * 1000
            s.RADr += r.PARi

           
    @prepare_states       
    def _set_variable_LAI(self,value):
        #print(f"Trying to update LAI with {value}")
        
        s = self.states
        oTLEAF=s.TLEAF        
        oLAI = max(0.1,s.LAI)
        nLAI=max(0.1,value)
        
        s.TLEAF=oTLEAF*nLAI/oLAI
        s.LAI=nLAI
        increments={"LAI":s.LAI - oLAI,"TLEAF":s.TLEAF-oTLEAF}
        
        
        return increments
