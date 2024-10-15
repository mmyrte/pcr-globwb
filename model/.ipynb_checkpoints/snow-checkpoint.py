import netCDF4 as nc
import pcraster as pcr
import numpy as np

import logging
logger = logging.getLogger(__name__)

import virtualOS as vos
from ncConverter import *
import glaciers as gl


def initializeSnow(self, iniItems):
    logger.info("Initialising Snow Extended...")
    #Initialize snow module.
    #First, a range of boolean parameters to specify how the snow module should work.
    
    #SnowTransport.
    if "snowTransport" in list(self.iniItemsLC.keys()):
        if self.iniItemsLC['snowTransport'] == "False": self.snowTransport = False
        else:
            self.snowTransport = self.iniItemsLC['snowTransport']
    else:
        self.snowTransport = False

    #Do we apply a precipitation/snowfall transition range?
    if "snowfallTransition" in list(self.iniItemsLC.keys()):
        if self.iniItemsLC['snowfallTransition'] == "True": self.snowfallTransition = True
        elif self.iniItemsLC['snowfallTransition'] == "False": self.snowfallTransition = False
    else:
        self.snowfallTransition = False

    #Do we use glaciers?
    if "glacierModule" in list(iniItems.landSurfaceOptions.keys()):
        if iniItems.landSurfaceOptions['glacierModule'] == "True": self.glacierModule = True
        elif iniItems.landSurfaceOptions['glacierModule'] == "False": self.glacierModule = False
    else:
        self.glacierModule = False

    #glacierType (static, deltaH, Immerzeel)  
    if "glacierType" in list(iniItems.landSurfaceOptions.keys()):
        self.glacierType = iniItems.landSurfaceOptions["glacierType"]
    else:
        self.glacierType = 'Static'

    #Exponential melt
    if "exponentialMelt" in list(self.iniItemsLC.keys()):
        if self.iniItemsLC['exponentialMelt'] == "True": self.exponentialMelt = True
        elif self.iniItemsLC['exponentialMelt'] == "False": self.exponentialMelt = False
    else:
        self.exponentialMelt = False

    #Albedo melt
    if "albedoMelt" in list(self.iniItemsLC.keys()):
        if self.iniItemsLC['albedoMelt'] == "True": self.albedoMelt = True
        elif self.iniItemsLC['albedoMelt'] == "False": self.albedoMelt = False
    else:
        self.albedoMelt = False

    #Seasonal melt
    if "seasonalMelt" in list(self.iniItemsLC.keys()):
        if self.iniItemsLC['seasonalMelt'] == "True": self.seasonalMelt = True
        elif self.iniItemsLC['seasonalMelt'] == "False": self.seasonalMelt = False
    else:
        self.seasonalMelt = False

    #Do we want to apply a snow fall correction?
    if "snowfallCorrection" in list(self.iniItemsLC.keys()):
        self.snowfallCorrection = float(self.iniItemsLC['snowfallCorrection'])
    else:
        self.snowfallCorrection=False

    #Precipitation driven melt
    if "precipitationMelt" in list(self.iniItemsLC.keys()):
        if self.iniItemsLC['precipitationMelt'] == "True": self.precipitationMelt = True
        elif self.iniItemsLC['precipitationMelt'] == "False": self.precipitationMelt = False
    else:
        self.precipitationMelt = False

    #Snow Module
    if self.albedoMelt:
        #Reading the parameters for albedo driven melt (Kraaijenbrink et al., 2020).
        logger.info("Initialising Albedo Melt...")
        self.DDFmin=self.degreeDayFactor
        
        input = self.iniItemsLC['degreeDayRange']
        vars(self)['degreeDayRange'] = vos.readPCRmapClone(input,self.cloneMap,
                                        self.tmpDir,self.inputDir)
        vars(self)['degreeDayRange'] = pcr.spatial(pcr.scalar(vars(self)['degreeDayRange']))
        
        self.DDFmax=self.DDFmin+self.degreeDayRange
        self.PDD=pcr.scalar(0)
        
    elif self.seasonalMelt:
        #Reading the parameters for seasonal melt (Slater and Clark, 2006).
        logger.info("Initialising Seasonal Melt...")
        self.DDFmin=self.degreeDayFactor
        input = self.iniItemsLC['degreeDayAmplitude']
        vars(self)['degreeDayAmplitude'] = vos.readPCRmapClone(input,self.cloneMap,
                                        self.tmpDir,self.inputDir)
        vars(self)['degreeDayAmplitude'] = pcr.spatial(pcr.scalar(vars(self)['degreeDayAmplitude']))

    else:
        #Rename degreeDayFactor to DDFmin for flexibility in the snow melt procedure later.
        self.DDFmin=self.degreeDayFactor
            
    #Read cell area; Cell area is also needed without snow transport
    self.cellArea = vos.readPCRmapClone(\
                  iniItems.routingOptions['cellAreaMap'],
                  self.cloneMap,self.tmpDir,self.inputDir)
    
    if self.snowTransport!=False:
        logger.info("Initialising Snow Transport...")
        #Read surface gradient
        self.gradient = vos.readPCRmapClone(iniItems.routingOptions[str('gradient')],\
                         self.cloneMap,self.tmpDir,self.inputDir)

        #Calculate distances within the gridcells.
        self.cellSizeInArcDeg = vos.getMapAttributes(self.cloneMap,"cellsize")  
        cellSizeInArcMin    =  self.cellSizeInArcDeg*60.
        self.verticalSizeInMeter =  cellSizeInArcMin*1852.  
        self.zonalDistance=self.cellArea/self.verticalSizeInMeter
        self.verticalSizeInMeter=float(self.verticalSizeInMeter)

        #Read LDD Map
        self.lddMap = vos.readPCRmapClone(iniItems.routingOptions['lddMap'],\
                                        self.cloneMap,self.tmpDir,self.inputDir,True)
        self.lddMap = pcr.lddrepair(pcr.ldd(self.lddMap))
        self.lddMap = pcr.lddrepair(self.lddMap)

        #Read DEM
        try: self.highResolutionDEM = vos.readPCRmapClone(\
                            iniItems.meteoDownscalingOptions['highResolutionDEM'],
                                   self.cloneMap,self.tmpDir,self.inputDir)
        except: self.highResolutionDEM = vos.readPCRmapClone(\
                                   iniItems.meteoOptions['highResolutionDEM'],
                                   self.cloneMap,self.tmpDir,self.inputDir)
        
        #Reading parameters for Frey and Holzmann.
        logger.info("Initialising Frey and Holzmann...")
        input = self.iniItemsLC['Hv']
        vars(self)['Hv'] = vos.readPCRmapClone(input,self.cloneMap,
                                        self.tmpDir,self.inputDir)
        vars(self)['Hv'] = pcr.spatial(pcr.scalar(vars(self)['Hv']))
        
        input = self.iniItemsLC['frho']
        vars(self)['frho'] = vos.readPCRmapClone(input,self.cloneMap,
                                        self.tmpDir,self.inputDir)
        vars(self)['frho'] = pcr.spatial(pcr.scalar(vars(self)['frho']))

        #Read surface slope and compute surface angles from this.
        slope=vos.readPCRmapClone(\
                        'global_30sec/landSurface/topography/merit_dem_processed/version_2021-02-XX/maps_covered_with_zero/tanslope_topography_parameters_30sec_february_2021_global_covered_with_zero.nc',
                               self.cloneMap,self.tmpDir,self.inputDir)
        logger.info('Compute new surface angles....')
        pcr.setglobaloption("degrees")
        self.angle=pcr.scalar(pcr.atan(slope))
        pcr.setglobaloption("radians")
        mn = pcr.cellvalue(pcr.mapminimum(self.angle),1)[0]
        mx = pcr.cellvalue(pcr.mapmaximum(self.angle),1)[0]
        logger.info("Mx: "+str(mx)+"; mn: "+str(mn))
    
        #Reverse DEM (needed to transport snow into multiple directions downslope).
        reverse_dem= self.highResolutionDEM*-1
        self.reverseLDD=pcr.lddcreate(reverse_dem,1e31,1e31,1e31,1e31)
        ones=pcr.scalar(1.0)
        self.downstreamCells=pcr.upstream(self.reverseLDD, ones)
        self.downstreamCells=pcr.ifthenelse(self.downstreamCells==0, 1.0, self.downstreamCells)
        #%%ADDED BY JOREN:STOP
    else:
     #For the water balance check, we need to define some empty stores if snowTransport is not used
        self.transportVolSnow=pcr.scalar(0.0)
        self.incomingVolSnow=pcr.scalar(0.0)

    #For the water balance check, we need to define some empty stores if glacierModule is not used:
    if self.glacierModule==False:
        self.glacierIce=pcr.scalar(0)
        self.glacierWater=pcr.scalar(0)

        self.glacierOutflow=pcr.scalar(0)
        self.netGlacierWaterToSoil=pcr.scalar(0)

        self.glacierAccumulation=pcr.scalar(0)
        self.rain2GlacierWater=pcr.scalar(0)
        self.snow2GlacierWater=pcr.scalar(0)
        self.capacity2GlacierWater=pcr.scalar(0)
        self.iceMelt=pcr.scalar(0)
        

def updateSnowFall(self, meteo,currTimeStep):
    #Apply a transition range between solid and liquid precipitation.
    logger.info('Snowfall transition range!')
    
    #Values from Magnusson et al., 2014
    self.Tpbase=1.0
    self.mp=1.24
    self.Tp=(meteo.temperature-self.Tpbase)/self.mp
    self.estimSnowfall=pcr.max(meteo.precipitation*1/(1+pcr.exp(self.Tp)), 0.0)

    if self.snowfallCorrection!=False:
        logger.info('-------Correcting snowfall on glaciers!!!--------------')
        self.estimSnowfall=pcr.ifthenelse(self.glacierIce>0, self.estimSnowfall*self.snowfallCorrection, self.estimSnowfall)

def snowMeltSlaterAndClark(self, meteo,currTimeStep):
    logger.info('Starting with the snow module!')
                    
    if self.debugWaterBalance:
        prevStates        = [self.snowCoverSWE,self.snowFreeWater, self.glacierIce, self.glacierWater]
        prevSnowCoverSWE  = self.snowCoverSWE
        prevSnowFreeWater = self.snowFreeWater
        prevGlacierModule = [self.glacierIce, self.glacierWater]
        prevGlacierWater  = self.glacierWater
        prevGlacierIce    = self.glacierIce
    
    #------SNOW-TRANSPORT-----------------------------------------------------------------
    #Snow transport can be perform wit
    self.transportVolSnow=pcr.scalar(0.0)
    self.incomingVolSnow=pcr.scalar(0.0)
    #What Kind Of Snow Transport Do We Want...;   
    #This is the default snow transport.
    if self.snowTransport=='FreyAndHolzmann_pcraster':
        simplifiedFreyAndHolzmann_pcraster(self, currTimeStep)
    else:
        logger.info('NO SnowSlide: who needs that anyway....?')
     #--------------------------------------------------------------------------------------   
        
    #------SNOW-MELT-----------------------------------------------------------------
    if self.albedoMelt==True:
        #Relate DDF to albedo of snow. (Kraaijenbrink et al., 2020)
        logger.info('Using Albedo Based Melt.....')
        self.PDD=pcr.ifthenelse(pcr.pcrand(self.snowfall<0.005, self.snowCoverSWE>0.0), self.PDD+pcr.max(meteo.temperature_max, 0), 0.001)
        self.albedo=0.713-0.112*pcr.log10(self.PDD)
        self.albedo=pcr.max(self.albedo, 0.4)
        self.albedo=pcr.min(self.albedo, 0.85)
        self.degreeDayFactor=(self.DDFmax-self.DDFmin)/(1/0.4-1/0.85)*(1/self.albedo-1/0.85)+self.DDFmin
    else:
        #Just use the given DDF.
        self.degreeDayFactor=self.DDFmin

    #Include melt due to precipitation.
    if self.precipitationMelt==True:
        logger.info('Using Precipitation Melt.....')
        Cp=4.18 #J/g/K
        Lf=333.55#J/g
        Pmelt=self.liquidPrecip*meteo.temperature*Cp/Lf
    else:
        Pmelt=pcr.scalar(0.0)

    
    if self.seasonalMelt==True:
        logger.info('Using Seasonal Melt.....')
        if self.exponentialMelt==False:
            #Calculate melt with seasonally varying DDFs, no exponential relationship to temperature.
            deltaSnowCover = \
                pcr.ifthenelse(meteo.temperature <= self.freezingT, \
                self.refreezingCoeff*self.snowFreeWater, \
               -pcr.min(self.snowCoverSWE, \
                        Pmelt+pcr.max(meteo.temperature - self.freezingT, 0.0) * pcr.max(self.degreeDayFactor+self.degreeDayAmplitude*np.sin((currTimeStep.doy-81)*2*np.pi/366), 0.0)))
        elif self.exponentialMelt==True:
            #Calculate melt with seasonally varying DDFs, this time with exponential relationship to temperature.
            Mm=0.5
            ratio=(meteo.temperature - self.freezingT)/Mm
            deltaSnowCover = -pcr.min(self.snowCoverSWE+self.snowfall, \
                        pcr.max(Pmelt+ Mm * (ratio + pcr.ln(1+pcr.exp(-1*ratio)))  * pcr.max(self.degreeDayFactor+self.degreeDayAmplitude*np.sin((currTimeStep.doy-81)*2*np.pi/366), 0.0), 0.0))
            #Should also melt new snowfall on glaciers!!!
    else:
        #Calculate melt with constant DDF.
        deltaSnowCover = \
            pcr.ifthenelse(meteo.temperature <= self.freezingT, \
            self.refreezingCoeff*self.snowFreeWater, \
           -pcr.min(self.snowCoverSWE, \
                    pcr.max(meteo.temperature - self.freezingT, 0.0) * self.degreeDayFactor+Pmelt))


    
    # update snowCoverSWE
    self.snowCoverSWE  = pcr.max(0.0, self.snowfall + deltaSnowCover + self.snowCoverSWE)                              
                                                                    # SC_L[TYPE] = max(0.0, SC_L[TYPE]+DSC[TYPE]+SNOW)

    # for reporting snow melt in m/day
    self.snowMelt = pcr.ifthenelse(deltaSnowCover < 0.0, deltaSnowCover * pcr.scalar(-1.0), pcr.scalar(0.0))
    
    #--------------------------------------------------------------------------------------
    
    

    #-----GLACIER-MODULE------------------------------------------------------------------
    if self.glacierModule==True:
        #Update static part of the glaciers: glacierMelt and accumulation.
        logger.info('Starting with the glaciers.....')
        self.deltaSnowCover=deltaSnowCover
        gl.updateStaticGlacier(self, meteo, currTimeStep)
        deltaSnowCover=self.deltaSnowCover
    
        self.netLqWaterToSoil = pcr.max(0., self.snowFreeWater - \
                 self.snowWaterHoldingCap * self.snowCoverSWE)
        
        # update snowFreeWater (after netLqWaterToSoil) 
        self.snowFreeWater    = pcr.max(0., self.snowFreeWater - \
                                            self.netLqWaterToSoil)      # SCF_L[TYPE] = max(0,SCF_L[TYPE]-Pn)
        
        self.netLqWaterToSoil=self.netLqWaterToSoil+self.netGlacierWaterToSoil
        
    # update snowFreeWater = liquid water stored above snowCoverSWE
    else:
        self.snowFreeWater = self.snowFreeWater - deltaSnowCover + \
                             self.liquidPrecip                          # SCF_L[TYPE] = SCF_L[TYPE]-DSC[TYPE]+PRP;                       
        # netLqWaterToSoil = net liquid transferred to soil
        self.netLqWaterToSoil = pcr.max(0., self.snowFreeWater - \
                 self.snowWaterHoldingCap * self.snowCoverSWE)          # Pn = max(0,SCF_L[TYPE]-CWH*SC_L[TYPE])
        
        # update snowFreeWater (after netLqWaterToSoil) 
        self.snowFreeWater    = pcr.max(0., self.snowFreeWater - \
                                            self.netLqWaterToSoil)      # SCF_L[TYPE] = max(0,SCF_L[TYPE]-Pn)
    #--------------------------------------------------------------------------------
    # evaporation from snowFreeWater (based on potBareSoilEvap)
    self.actSnowFreeWaterEvap = pcr.min(self.snowFreeWater, \
                                        self.potBareSoilEvap)       # ES_a[TYPE] = min(SCF_L[TYPE],ES_p[TYPE])
                                   
    # update snowFreeWater and potBareSoilEvap
    self.snowFreeWater = pcr.max(0.0, \
                         self.snowFreeWater - self.actSnowFreeWaterEvap)  
                                                                    # SCF_L[TYPE]= SCF_L[TYPE]-ES_a[TYPE]
    self.potBareSoilEvap = pcr.max(0, \
                       self.potBareSoilEvap - self.actSnowFreeWaterEvap) 
                                                                    # ES_p[TYPE]= max(0,ES_p[TYPE]-ES_a[TYPE])

    # update actual evaporation (after evaporation from snowFreeWater) 
    self.actualET += self.actSnowFreeWaterEvap                      # EACT_L[TYPE]= EACT_L[TYPE]+ES_a[TYPE];
    


    #Perform water balance check.
    if self.debugWaterBalance:
        if self.glacierModule==True:
            logger.info('WBCHECK including GlacierModule')
        else:
            logger.info('WBCHECK NO GlacierModule')
            self.glacierIce=pcr.scalar(0)
            self.glacierWater=pcr.scalar(0)

            self.glacierOutflow=pcr.scalar(0)
            self.netGlacierWaterToSoil=pcr.scalar(0)

            self.glacierAccumulation=pcr.scalar(0)
            self.rain2GlacierWater=pcr.scalar(0)
            self.snow2GlacierWater=pcr.scalar(0)
            self.capacity2GlacierWater=pcr.scalar(0)
            self.iceMelt=pcr.scalar(0)
                    
        if self.snowTransport=='Complete' or self.snowTransport=='Basic' or self.snowTransport=='FreyAndHolzmann' or self.snowTransport=='FreyAndHolzmann_pcraster':
            logger.info('WBCHECK including SnowTransport: '+str(self.snowTransport))
        else:
            logger.info('WBCHECK including SnowTransport: '+str(self.snowTransport))
            self.incomingVolSnow=pcr.scalar(0)
            self.transportVolSnow=pcr.scalar(0)
            
        
        vos.waterBalanceCheck([self.snowfall, self.liquidPrecip, self.incomingVolSnow/self.cellArea],
                              [self.netLqWaterToSoil,\
                               self.actSnowFreeWaterEvap, self.transportVolSnow/self.cellArea, self.glacierOutflow],
                               prevStates,\
                              [self.snowCoverSWE, self.snowFreeWater, self.glacierIce, self.glacierWater],\
                              'snow module local',\
                               True,\
                               currTimeStep.fulldate,threshold=1e-4)

        vos.waterBalanceCheck([self.snowfall, deltaSnowCover, self.incomingVolSnow/self.cellArea],\
                              [self.transportVolSnow/self.cellArea, self.glacierAccumulation],\
                              [prevSnowCoverSWE],\
                              [self.snowCoverSWE],\
                              'snowCoverSWE',\
                               True,\
                               currTimeStep.fulldate,threshold=1e-4)

        vos.waterBalanceCheck([self.liquidPrecip-self.rain2GlacierWater, -1*deltaSnowCover],
                              [self.actSnowFreeWaterEvap, self.netLqWaterToSoil-self.netGlacierWaterToSoil, self.snow2GlacierWater, self.capacity2GlacierWater],
                              [prevSnowFreeWater],\
                              [self.snowFreeWater],\
                              'snowFreeWater',\
                               True,\
                               currTimeStep.fulldate,threshold=1e-4)
        

        vos.waterBalanceCheck([self.rain2GlacierWater, self.snow2GlacierWater, self.glacierAccumulation, self.capacity2GlacierWater],
                              [self.glacierOutflow, self.netGlacierWaterToSoil],
                              prevGlacierModule,\
                              [self.glacierIce, self.glacierWater],\
                              'glacierModule local',\
                               True,\
                               currTimeStep.fulldate,threshold=1e-4)
            
            
            
        vos.waterBalanceCheck([self.rain2GlacierWater, self.snow2GlacierWater, self.iceMelt, self.capacity2GlacierWater],
                              [self.glacierOutflow, self.netGlacierWaterToSoil],
                              [prevGlacierWater],\
                              [self.glacierWater],\
                              'glacierWater',\
                               True,\
                               currTimeStep.fulldate,threshold=1e-4)
        
        
        vos.waterBalanceCheck([self.glacierAccumulation],
                              [self.iceMelt],
                              [prevGlacierIce],\
                              [self.glacierIce],\
                              'glacierIce local',\
                               True,\
                               currTimeStep.fulldate,threshold=1e-4)




def simplifiedFreyAndHolzmann_pcraster(self, currTimeStep):
    logger.info('Starting with Frey and Holzmann PCRASTER: let\'s see how far we get....')
    
    #TRANSPORT SNOW
    #Check where snow exceeds the threshold.
    exceedingSnow=pcr.max(self.snowCoverSWE-self.Hv, 0.0)
    #if self.glacierModule==True:
        #if currTimeStep.doy==275:
        #if (currTimeStep.doy>=275) & (currTimeStep.doy<=300):
        #    logger.info('Today with transport on Glaciers!!....')
        #    print('Today with transport on Glaciers!!....')
        #else:
        #    logger.info('No Transport on Glaciers!!....')
        #    exceedingSnow=pcr.ifthenelse(self.glacierized>0, 0.0, exceedingSnow)
    
    #Convert everything to volumes
    transportVolSnow=pcr.max(exceedingSnow*self.cellArea, 0.0)
    
    #Calculate fraction that needs to be transported (based on surface slope)
    self.transportVolSnow=transportVolSnow*self.angle/90*self.frho
    
    #Divide by number of downstream cells (downstream copies the value of the downstream cell, thus this step is needed for water balance)        
    fractionTransport=self.transportVolSnow/self.downstreamCells
    #Transport the snow to downstream cells (with reverse LDD)
    self.incomingVolSnow = pcr.downstream(self.reverseLDD, fractionTransport)
    
    #TRANSPORT SNOWFREEWATER
    self.transport_water=False
    if self.transport_water==True:
        frac_of_snow=(exceedingSnow*self.angle/90*self.frho)/self.snowCoverSWE
        self.transportFreeWater=pcr.max(self.cellArea*self.snowFreeWater*frac_of_snow, 0.0)
        fractionWater=self.transportFreeWater/self.downstreamCells
        #Transport the snow free water to downstream cells (with reverse LDD)
        self.incomingFreeWater = pcr.downstream(self.reverseLDD, fractionWater)
        self.snowFreeWater = self.snowFreeWater-self.transportFreeWater/self.cellArea+self.incomingFreeWater/self.cellArea
        
    #Compute new snow cover
    self.snowCoverSWE = self.snowCoverSWE-self.transportVolSnow/self.cellArea+self.incomingVolSnow/self.cellArea
    
    logger.info('Finished with Frey and Holzmann PCRASTER! Impressive..')
            



















