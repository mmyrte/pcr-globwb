import netCDF4 as nc
import pcraster as pcr
import numpy as np

import logging
logger = logging.getLogger(__name__)

import virtualOS as vos
from ncConverter import *


#Still needs to be added:
# - SnowFreeWater exceeding capacity added to the glacier water
# - Glacier water needs to have a max capacity?
# - Glacier water should refreeze?
# - Glacier volume?

def updateStaticGlacier(self, meteo, currTimeStep):
    #To Do: make sure it only happens if there is a glacier.
    
    #Empirical constants
    #For melting:
    MRG=1.52 #Melt Rate of Glacier over Melt Rate of Snow; Stahl et al., 2008

    #For outflow:
    Kmin=0.01 #Minimal outflow fraction
    Krange=0.03 #Range
    Ag=2 #Dependency on snow cover
    
    #For accummulation
    AccEfficiency=0.001#Fraction of snow turned to ice
    
    #Glacier capacity
    glacierWaterHoldingCap=self.snowWaterHoldingCap
    
    
    #Calculate Glacier Accummulation
    self.glacierAccummulation=AccEfficiency*self.glacierized*self.snowCoverSWE
    self.glacierIce+=self.glacierAccummulation
    self.snowCoverSWE-=self.glacierAccummulation
    
        
    #Ice Melt
    #To Do: make melt different from snow; Add rain?
    self.iceMelt = \
        pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingT),(self.snowCoverSWE>0.0)), \
        0.0, \
       -pcr.min(self.glacierIce, \
                self.glacierized*pcr.max(meteo.temperature - self.freezingT, 0.0) * MRG * (self.degreeDayFactor+self.degreeDayAmplitude*np.sin((currTimeStep.doy-81)*2*np.pi/366))))
    
    self.glacierIce-=self.iceMelt
    
    
    #Update Glacier Water
    self.snow2GlacierWater= - pcr.min(0.0, self.deltaSnowCover) * self.glacierized #Convert to positive when added to glacier, remove refreezing
      
    self.rain2GlacierWater=pcr.ifthenelse(pcr.pcrand(self.glacierized>0.0, self.snowCoverSWE==0), self.liquidPrecip, 0.0)

    self.snowFreeWater = self.snowFreeWater \
                                 - self.deltaSnowCover - self.snow2GlacierWater \
                                 +  self.liquidPrecip -  self.rain2GlacierWater
    
    self.capacity2GlacierWater=pcr.ifthenelse(self.glacierized>0.0, pcr.max(0., self.snowFreeWater - self.snowWaterHoldingCap * self.snowCoverSWE), 0.0)
   
    # update snowFreeWater (after capacity2GlacierWater) 
    self.snowFreeWater    = pcr.max(0., self.snowFreeWater - \
                                            self.capacity2GlacierWater)

    self.glacierWater=self.glacierWater+self.rain2GlacierWater+self.snow2GlacierWater+self.capacity2GlacierWater+self.iceMelt
    
    #TO DO:
    #Add overcapacity to soil? Or to outflow?
    #Add own capicity for glacier
    #Add refreezing?
    
    self.netGlacierWaterToSoil = pcr.max(0., self.glacierWater - \
                     glacierWaterHoldingCap * self.glacierIce)
    
    self.glacierWater=self.glacierWater-self.netGlacierWaterToSoil
    
    
    
    #Calculate Outflow and add it to the snowmelt.
    #TODO: get values from literature: Values are in Van Tiel et al., 2018
    self.glacierOutflow=self.glacierWater*(Kmin+Krange*pcr.exp(-Ag*self.snowCoverSWE))*self.glacierized
    
    
    self.glacierWater=self.glacierWater-self.glacierOutflow
    #self.snowFreeWater=self.snowFreeWater+self.glacierOutflow
    
    self.directRunoff=self.directRunoff+self.glacierOutflow