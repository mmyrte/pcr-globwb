import netCDF4 as nc
import pcraster as pcr
import numpy as np

import logging
logger = logging.getLogger(__name__)

import virtualOS as vos
from ncConverter import *


#Still needs to be added:
# - SnowFreeWater exceeding capacity added to the glacier water
# - Make sure they can't be negative!!!!!!!!
# - Glacier water needs to have a max capacity?
# - Glacier water should refreeze?
# - Glacier volume?

def initializeGlacier(self, iniItems):
    logger.info("Initialising Glaciers...")
    self.glacierIce = vos.readPCRmapClone(\
                      iniItems.landSurfaceOptions['glaciers'],
                      self.cloneMap,self.tmpDir,self.inputDir)
    self.glacierIce=pcr.ifthenelse(self.glacierIce>0.0, self.glacierIce, 0.0)
    self.glacierized=pcr.scalar(self.glacierIce>0.0)
    self.glacierWater=pcr.ifthenelse(self.glacierized>0, 0.01*self.glacierIce, 0)
    self.glacierInitial=self.glacierIce
    self.gradient = vos.readPCRmapClone(iniItems.routingOptions[str('gradient')],\
                             self.cloneMap,self.tmpDir,self.inputDir)
    
    
    #Values from Immerzeel et al., 2012
    #tau_0=80000 #Equilibrium shear stress [Nm-2]
    #rho=916.7 #Ice density [kgm-3]
    #R=1e-9 #Material roughness coefficient [Nm-2s1/3]
    #nu=0.1 #Bedrock roughness [-]
    #g=9.81 #gravitational acceleration [ms-2]
    #n=3 #creep constant from Glenns law

    #H_balance=tau_0/(rho*g*pcr.sin(pcr.scalar(pcr.atan(self.gradient))))
    #self.glacierIce=pcr.ifthenelse(self.glacierIce>H_balance, H_balance+1, self.glacierIce)
    #logger.info("Reshaping Glaciers...")
    
def updateStaticGlacier(self, meteo, currTimeStep):
    #To Do: make sure it only happens if there is a glacier.
    
    #Empirical constants
    #For melting:
    MRG=1 #Melt Rate of Glacier over Melt Rate of Snow; Stahl et al., 2008
    self.degreeDayFactorGlacier=0.009
    
    #For outflow:
    #Derived from Van Tiel et al., 2018
    Kmin=0.2 #Minimal outflow fraction
    Krange=0.5 #Range
    Ag=1.25 #Dependency on snow cover
    
    #For accummulation
    AccEfficiency=0.001#Fraction of snow turned to ice
    
    #Glacier capacity
    glacierWaterHoldingCap=self.snowWaterHoldingCap
    acclim=0.625 #minimum snow thickness needed for accummulation [m]
    
    #Calculate Glacier Accummulation
    self.glacierAccummulation=pcr.ifthenelse(self.snowCoverSWE>acclim, AccEfficiency*self.glacierized*self.snowCoverSWE, 0)
    self.glacierIce+=self.glacierAccummulation
    self.snowCoverSWE-=self.glacierAccummulation
    
        
    #Ice Melt
    #To Do: make melt different from snow; Add rain?
    if self.seasonalMelt:
        self.iceMelt = \
            pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingT),(self.snowCoverSWE>0.0)), \
            0.0, \
           pcr.min(self.glacierIce, \
                    self.glacierized*pcr.max(meteo.temperature - self.freezingT, 0.0) * MRG * (self.degreeDayFactor+self.degreeDayAmplitude*np.sin((currTimeStep.doy-81)*2*np.pi/366))))
    else:
        self.iceMelt = \
            pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingT),(self.snowCoverSWE>0.0)), \
            0.0, \
           pcr.min(self.glacierIce, \
                    self.glacierized*pcr.max(meteo.temperature - self.freezingT, 0.0) * MRG * (self.degreeDayFactorGlacier)))
        
    
    self.glacierIce-=self.iceMelt
    
    
    #Update Glacier Water
    self.snow2GlacierWater= - pcr.min(0.0, self.deltaSnowCover) * self.glacierized #Convert to positive when added to glacier, remove refreezing
      
    self.rain2GlacierWater=pcr.ifthenelse(pcr.pcrand(self.glacierized>0.0, self.snowCoverSWE==0), self.liquidPrecip, 0.0)

    self.snowFreeWater = pcr.max(0.0, self.snowFreeWater \
                                 - self.deltaSnowCover - self.snow2GlacierWater \
                                 +  self.liquidPrecip -  self.rain2GlacierWater)
    
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
    
    self.glacierWater=pcr.max(0, self.glacierWater-self.netGlacierWaterToSoil)
    
    
    
    #Calculate Outflow and add it to the snowmelt.
    #TODO: get values from literature: Values are in Van Tiel et al., 2018
    self.glacierOutflow=pcr.min(self.glacierWater, self.glacierWater*(Kmin+Krange*pcr.exp(-Ag*self.snowCoverSWE))*self.glacierized)
    
    
    self.glacierWater=pcr.max(0, self.glacierWater-self.glacierOutflow)
    #self.snowFreeWater=self.snowFreeWater+self.glacierOutflow
    

    

def glacierSlideImmerzeel(self, currTimeStep):
    logger.info('Starting with GlacierSlideImmerzeel: let\'s see how far we get....')
    #critSWE=self.maxSWE(angle)

    #Values from Immerzeel et al., 2012
    tau_0=80000 #Equilibrium shear stress [Nm-2]
    rho=916.7 #Ice density [kgm-3]
    R=1e9 #Material roughness coefficient [Nm-2s1/3]
    nu=0.1 #Bedrock roughness [-]
    g=9.81 #gravitational acceleration [ms-2]
    n=3 #creep constant from Glenns law

    flowlim=5 #minimum ice thickness needed for flow [m]

    #Copied from routing; as a first guess
    # channelLength = approximation of channel length (unit: m)
    # This is approximated by cell diagonal.                     
    self.cellLengthFD  = ((self.cellArea/self.verticalSizeInMeter)**(2)+\
                                        (self.verticalSizeInMeter)**(2))**(0.5) 
    self.channelLength = self.cellLengthFD


    u_two_fourth=pcr.max(0,(rho*g*self.glacierIce*pcr.sin(pcr.scalar(pcr.atan(self.gradient))) - tau_0)/(nu**2*R)) #[(m/s)**1/2]
    u=u_two_fourth**((n+1)/2) #[m/s]

    self.glacierVelocity=pcr.max(0,(u*3600*24)) #[m/day]
    a,b,c =vos.getMinMaxMean(self.glacierVelocity)
    msg = "Glacier Velocities initialized: Min %f Max %f Mean %f" %(a,b,c)
    logger.info(msg)
    
    #self.glacierVelocity=pcr.min(self.glacierVelocity, 3.5)
    a,b,c =vos.getMinMaxMean(self.glacierVelocity)
    msg = "Glacier Velocities corrected: Min %f Max %f Mean %f" %(a,b,c)
    logger.info(msg)
                                 
    #U=self.glacierVelocity/self.channelLength #[cells/day]

    # updating (after routing)
    self.glacierIceVolume=self.glacierIce*self.cellArea
    self.glacierOutgoing = pcr.ifthenelse(self.glacierIce>flowlim,self.glacierVelocity*(self.glacierIceVolume/self.channelLength), 0)
    self.glacierOutgoing = pcr.ifthenelse(self.glacierOutgoing>self.glacierIceVolume, self.glacierIceVolume, self.glacierOutgoing)
    self.glacierIncoming = pcr.upstream(self.lddMap,\
                          self.glacierOutgoing)
    self.glacierIceVolume=pcr.max(0, self.glacierIceVolume+self.glacierIncoming-self.glacierOutgoing)
    #self.glacierDelta=pcr.ifthenelse(self.glacierIceVolume+self.glacierIncoming-self.glacierOutgoing>0,0, self.glacierIceVolume+self.glacierIncoming-self.glacierOutgoing)
    #self.glacierOutgoing-=self.glacierDelta
    self.glacierIce=self.glacierIceVolume/self.cellArea
    
    self.glacierized=pcr.scalar(self.glacierIce>0.0)
    #self.glacierIce=self.glacierIce*self.glacierized

    logger.info('Finished with GlacierSlideImmerzeel: Impressive....')