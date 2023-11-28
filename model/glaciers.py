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
    self.glacierIce_ini=self.glacierIce
    
    self.glacierized=pcr.scalar(self.glacierIce>0.0)
    self.glacierWater=pcr.ifthenelse(self.glacierized>0, 0.01*self.glacierIce, 0)
    
    self.gradient = vos.readPCRmapClone(iniItems.routingOptions[str('gradient')],\
                             self.cloneMap,self.tmpDir,self.inputDir)
        
    self.glacierNumber = vos.readPCRmapClone(\
                      iniItems.landSurfaceOptions['glacierNumber'],
                      self.cloneMap,self.tmpDir,self.inputDir)
    
    #self.glacierNumber = pcr.readmap('/home/gjanzing/data/glaciers_schweiz_correct_new.map') #xr.open_dataset('/home/gjanzing/data/glaciers_schweiz_correct_new.nc')
    self.glacierNumber = pcr.scalar(self.glacierNumber)
    
    
    #----Read Variables for updating glaciers:--------------------------------------------------------
    glacierParams      = ['degreeDayFactorGlacier','acclim','Kmin','Krange','Ag'] #'AccEfficiency',
    
    for var in glacierParams:
        input = iniItems.landSurfaceOptions[str(var)]
        vars(self)[var] = vos.readPCRmapClone(input,self.cloneMap,
                                        self.tmpDir,self.inputDir)
        vars(self)[var] = pcr.spatial(pcr.scalar(vars(self)[var]))
    
    for coverType in self.coverTypes:
        logger.info("Dividing glacier parameters over: "+str(coverType))
        for var in glacierParams:
            setattr(self.landCoverObj[coverType], var, vars(self)[var])
        
    
    if self.glacierType=='delta_H':
        #Keep track of the amount of melt, Accumulation and the glacier ice at the start of the year.
        self.yearlyIceMelt=pcr.scalar(0.0)
        self.yearlyGlacierAcc=pcr.scalar(0.0)
        self.glacierIceYS=self.glacierIce
        
        if "ini_dH" in list(self.iniItems.landSurfaceOptions.keys()):
            if self.iniItems.landSurfaceOptions['ini_dH'] == "True": self.ini_dH = True
            elif self.iniItems.landSurfaceOptions['ini_dH'] == "False": self.ini_dH = False
        else:
            self.ini_dH = False
        
        if self.ini_dH==True:
            gl.initializeDeltaH(self)
    
    if self.glacierType=='Immerzeel':
        
        self.tau_0=364000 #Equilibrium shear stress [Nm-2]
        self.R=5.8e10 #Material roughness coefficient [Nm-2s1/3]
        self.flowlim=5 #[m]; Flow limit. If ice is thinner, there is no flow
        
        self.lddMap = vos.readPCRmapClone(iniItems.routingOptions['lddMap'],\
                                            self.cloneMap,self.tmpDir,self.inputDir,True)
        self.lddMap = pcr.lddrepair(pcr.ldd(self.lddMap))
        self.lddMap = pcr.lddrepair(self.lddMap)
        
        
        self.cellSizeInArcDeg = vos.getMapAttributes(self.cloneMap,"cellsize")  
        cellSizeInArcMin    =  self.cellSizeInArcDeg*60.
        self.verticalSizeInMeter =  cellSizeInArcMin*1852.  

        self.zonalDistance=self.cellArea/self.verticalSizeInMeter
        self.verticalSizeInMeter=float(self.verticalSizeInMeter)

        #Copied from routing; as a first guess
        # channelLength = approximation of channel length (unit: m)
        # This is approximated by cell diagonal.                     
        self.cellLengthFD  = ((self.cellArea/self.verticalSizeInMeter)**(2)+\
                                            (self.verticalSizeInMeter)**(2))**(0.5) 
        self.channelLength = self.cellLengthFD
    
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
    
def updateStaticGlacier_old(self, meteo, currTimeStep):
    #To Do: make sure it only happens if there is a glacier.
    
    #Empirical constants
    #For melting:
    #MRG=1 #Melt Rate of Glacier over Melt Rate of Snow; Stahl et al., 2008
    #self.degreeDayFactorGlacier=0.0055
    
    #For outflow:
    #Derived from Van Tiel et al., 2018
    #Kmin=0.2 #Minimal outflow fraction
    #Krange=0.5 #Range
    #Ag=1.25 #Dependency on snow cover
    
    #For Accumulation
    #AccEfficiency=0.005#Fraction of snow turned to ice
    
    #Glacier capacity
    glacierWaterHoldingCap=self.snowWaterHoldingCap
    #acclim=0 #minimum snow thickness needed for Accumulation [m]
    
    #Calculate Glacier Accumulation
    self.glacierAccumulation=pcr.ifthenelse(self.snowCoverSWE>self.acclim, self.AccEfficiency*self.glacierized*self.snowCoverSWE, 0)
    self.glacierIce+=self.glacierAccumulation
    self.snowCoverSWE-=self.glacierAccumulation
    
        
    #Ice Melt
    #To Do: make melt different from snow; Add rain?
    #if self.seasonalMelt:
    #    self.iceMelt = \
    #        pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingT),(self.snowCoverSWE>0.0)), \
    #        0.0, \
    #       pcr.min(self.glacierIce, \
    #                self.glacierized*pcr.max(meteo.temperature - self.freezingT, 0.0)*(self.degreeDayFactor+self.degreeDayAmplitude*np.sin((currTimeStep.doy-81)*2*np.pi/366))))
    #else:
    self.iceMelt = \
        pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingT),(self.snowCoverSWE>0.0)), \
        0.0, \
       pcr.min(self.glacierIce, \
                self.glacierized*pcr.max(meteo.temperature - self.freezingT, 0.0) * (self.degreeDayFactorGlacier)))
        
    
    self.glacierIce-=self.iceMelt
    
    if self.glacierType=='delta_H':
        self.yearlyIceMelt+=self.iceMelt
        self.yearlyGlacierAcc+=self.glacierAccumulation
            
        
        
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
    self.glacierOutflow=pcr.min(self.glacierWater, self.glacierWater*(self.Kmin+self.Krange*pcr.exp(-1*self.Ag*self.snowCoverSWE))*self.glacierized)
    
    
    self.glacierWater=pcr.max(0, self.glacierWater-self.glacierOutflow)
    #self.snowFreeWater=self.snowFreeWater+self.glacierOutflow
    
def updateStaticGlacier(self, meteo, currTimeStep):
    #To Do: make sure it only happens if there is a glacier.
    
    #Empirical constants
    #For melting:
    #MRG=1 #Melt Rate of Glacier over Melt Rate of Snow; Stahl et al., 2008
    #self.degreeDayFactorGlacier=0.0055
    
    #For outflow:
    #Derived from Van Tiel et al., 2018
    #Kmin=0.2 #Minimal outflow fraction
    #Krange=0.5 #Range
    #Ag=1.25 #Dependency on snow cover
    
    #For Accumulation
    #AccEfficiency=0.005#Fraction of snow turned to ice
    
    #Glacier capacity
    glacierWaterHoldingCap=self.snowWaterHoldingCap
    #acclim=0 #minimum snow thickness needed for Accumulation [m]
    
    #Calculate Glacier Accumulation
    if currTimeStep.doy==244:
        #self.glacierAccumulation=pcr.ifthenelse(self.snowCoverSWE>self.acclim, self.AccEfficiency*self.glacierized*self.snowCoverSWE, 0)
        self.glacierAccumulation=pcr.ifthenelse(self.snowCoverSWE>self.acclim, self.snowCoverSWE-self.acclim, 0)
        self.glacierIce+=self.glacierAccumulation
        self.snowCoverSWE-=self.glacierAccumulation
    else:
        self.glacierAccumulation=pcr.scalar(0.0)
    
        
    #Ice Melt
    #To Do: make melt different from snow; Add rain?
    #if self.seasonalMelt:
    #    self.iceMelt = \
    #        pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingT),(self.snowCoverSWE>0.0)), \
    #        0.0, \
    #       pcr.min(self.glacierIce, \
    #                self.glacierized*pcr.max(meteo.temperature - self.freezingT, 0.0)*(self.degreeDayFactor+self.degreeDayAmplitude*np.sin((currTimeStep.doy-81)*2*np.pi/366))))
    #else:
    
    self.iceMelt = \
        pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingT),(self.snowCoverSWE>0.0)), \
        0.0, \
       pcr.min(self.glacierIce, \
                self.glacierized*pcr.max(meteo.temperature - self.freezingT, 0.0) * (self.degreeDayFactorGlacier)))
        
    
    self.glacierIce-=self.iceMelt
    
    if self.glacierType=='delta_H':
        self.yearlyIceMelt+=self.iceMelt
        self.yearlyGlacierAcc+=self.glacierAccumulation
            
        
        
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
    self.glacierOutflow=pcr.min(self.glacierWater, self.glacierWater*(self.Kmin+self.Krange*pcr.exp(-1*self.Ag*self.snowCoverSWE))*self.glacierized)
    
    
    self.glacierWater=pcr.max(0, self.glacierWater-self.glacierOutflow)
    #self.snowFreeWater=self.snowFreeWater+self.glacierOutflow
    

def glacierSlideImmerzeel(self, currTimeStep):
    logger.info('Starting with GlacierSlideImmerzeel: let\'s see how far we get....')
    #critSWE=self.maxSWE(angle)

    #Values from Immerzeel et al., 2012
    tau_0=self.tau_0 #Equilibrium shear stress [Nm-2]
    rho=916.7 #Ice density [kgm-3]
    R=self.R #Material roughness coefficient [Nm-2s1/3]
    nu=0.1 #Bedrock roughness [-]
    g=9.81 #gravitational acceleration [ms-2]
    n=3 #creep constant from Glenns law

    flowlim=self.flowlim #minimum ice thickness needed for flow [m]

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
    

def updateDeltaH(self, currTimeStep):
    logger.info('Updating delta H.......')
    mn = pcr.cellvalue(pcr.mapminimum(self.glacierNumber),1)[0]
    mx = pcr.cellvalue(pcr.mapmaximum(self.glacierNumber),1)[0]

    glacierIcenew=pcr.scalar(0.0)
    extraMelt=pcr.scalar(0.0)
    
    self.glacierMB=self.yearlyGlacierAcc-self.yearlyIceMelt

    for i in range(int(mn), int(mx)+1):
    #for i in range(int(mn), int(mn)+30):
        glacierShape=self.glacierNumber==i
        volume_ini=pcr.ifthenelse(glacierShape, self.glacierIce_ini, 0.0)
        volume_ini=pcr.ifthenelse(volume_ini>0.0,volume_ini,0.0)
        glacierShape=volume_ini>0
        valid=pcr.cellvalue(pcr.maptotal(pcr.scalar(volume_ini>0)),1)[0]

        if valid==0:
            continue

        M=pcr.cellvalue(pcr.maptotal(volume_ini),1)[0]

        volume=pcr.ifthenelse(glacierShape, self.glacierIceYS, 0.0)
        volume=pcr.ifthenelse(volume>0.0,volume,0.0)
        Mpresent=pcr.cellvalue(pcr.maptotal(volume),1)[0]

        indivMB=pcr.ifthenelse(glacierShape, self.glacierMB, 0.0)
        deltaM=pcr.cellvalue(pcr.maptotal(indivMB),1)[0]

        frac_min=np.round((deltaM+Mpresent)/M,2)
        frac=1-frac_min
        if frac>1:
            frac=1

        #We don't want to have glacier growing when they are actually loosing mass due to round-off errors.
        #Therefore, in case a glacier does not move to another percentage point, we just take the old shape
        old_frac=1-np.round((Mpresent)/M,2)
        if int(frac*100)==int(old_frac*100):
            indivNew=self.glacierIceYS
        elif frac>0:
            indivNew=pcr.readmap("/hyclimm/gjanzing/data/glacierShape_minus"+str(int(frac*100)).zfill(2)+"M.map")
        elif frac<0:
            #print('Sir, we have caught a grower!!')
            indivNew=pcr.readmap("/hyclimm/gjanzing/data/glacierShape_minus00M.map")
        
        #FIX: cover NaN values
        indivNew=pcr.cover(indivNew, 0.0)
        
        newShape=pcr.ifthenelse(glacierShape, indivNew, 0.0)
        glacierIcenew=glacierIcenew+newShape

        Mnew=pcr.cellvalue(pcr.maptotal(newShape),1)[0]

        #To Do: what to do with extra volume??
        #--> add to of subtract from glacier Ice
        # else groundwater?

        #If difference between previous shape and new shape is larger than delta M (glacierDiff is positive),
        #this means that there is more melt than is observed.
        #Extra melt has to be added to the glacier.
        #And Vice Versa
        glacierGap=(Mpresent-Mnew)+deltaM

        #print('Mnew: '+str(Mnew))
        #print('actual delta: '+str(Mnew-Mpresent))

        #Only divide too much or too little melt over the glacier.
        #NOTE: also when a glacier grows more than the initial shape, it will be divided over the glacier as glacier ice!!

        if (frac!=1) and (abs(glacierGap)>1e-3):
            #We only divide stuff over places where the thickness is larger than 1, in order to prevent cells where the glacierDiff is larger than the ice thickness.
            area=pcr.cellvalue(pcr.maptotal(pcr.scalar(newShape>1)),1)[0]

            #If all cells are smaller than 1m, than area=1.
            if area==0:
                area=1

            #glacierDiff=pcr.ifthenelse(glacierShape, pcr.scalar(glacierGap/area), 0.0)
            glacierDiff=pcr.ifthenelse(newShape>1, pcr.scalar(glacierGap/area), 0.0)
            
            #Glaciers can not grow beyond their original extend. In case this happens, it is added to snow. 
            #This can then lead to snow towers, which could be removed by snow redistribution.
            #if (frac<0) or ((frac==0) and (glacierGap>0)):
            #    self.snowCoverSWE+=glacierDiff
                
            #    #This will be overwritten anyway, right?
            #    self.glacierAccuCorrection+=glacierDiff
            #else:
            glacierIcenew=glacierIcenew+glacierDiff
        
        
        #Should not happen: melt is restricted to the glacier volume, right?
        #Should this be though??
        extraMelt+=pcr.min(glacierIcenew, 0.0)
        glacierIcenew=pcr.max(glacierIcenew, 0.0)
        
        #Implement a good water balance check:..... but how?
        #assert pcr.cellvalue(pcr.mapminimum(pcr.scalar(extraMelt)),1)[0] >= -0.01
    
    em_error = pcr.cellvalue(pcr.mapminimum(extraMelt),1)[0]
    
    a,b,c =vos.getMinMaxMean(extraMelt)
    msg = "Extra melt: Min %f Max %f Mean %f" %(a,b,c)
    print(msg)

    if em_error<-1e-4:
        logger.info('ERROR: EXTRA MELT: '+str(em_error))

    #Update the new GlacierIce!
    self.glacierIce=glacierIcenew 
    #FIX: cover NaN values
    self.glacierIce=pcr.cover(self.glacierIce, 0.0)
    self.glacierized=pcr.scalar(self.glacierIce>0.0)

    #Re-initialize the storages over time.
    self.glacierIceYS=self.glacierIce
    self.yearlyIceMelt=pcr.scalar(0.0)
    self.yearlyGlacierAcc=pcr.scalar(0.0)
    logger.info('Finished with delta H!!.......')

def initializeDeltaH(self):
    logger.info('Initializing Delta H parameterization....')

    mn = pcr.cellvalue(pcr.mapminimum(self.glacierNumber),1)[0]
    mx = pcr.cellvalue(pcr.mapmaximum(self.glacierNumber),1)[0]


    #fs_list=np.arange(0.01, 1, 0.01)
    fs_list=np.arange(0, 101, 1)
    #fs_list=np.arange(0, 1+1e-5, 0.1)
    for frac in fs_list:
        print(frac)
        fs_map=pcr.scalar(0.0)
        logger.info("Creating: /hyclimm/gjanzing/data/glacierShape_minus"+str(int(frac)).zfill(2)+"M.map")
        pcr.report(fs_map,"/hyclimm/gjanzing/data/glacierShape_minus"+str(int(frac)).zfill(2)+"M.map")

    #mn=16823
    #print(mx)
    for i in range(int(mn), int(mx)+1):

        glacierShape=self.glacierNumber==i
        volume=pcr.ifthenelse(glacierShape, self.glacierIce, 0.0)
        volume=pcr.ifthenelse(volume>0.0,volume,0.0)
        glacierShape=volume>0
        valid=pcr.cellvalue(pcr.maptotal(pcr.scalar(volume>0)),1)[0]

        if valid==0:
            #When there are no cells, continue
            continue

        M=pcr.cellvalue(pcr.maptotal(volume),1)[0]

        glacierElev=pcr.ifthenelse(glacierShape, self.highResolutionDEM, np.nan)
        mn_elev = pcr.cellvalue(pcr.mapminimum(glacierElev),1)[0]
        mx_elev = pcr.cellvalue(pcr.mapmaximum(glacierElev),1)[0]

        elevs=np.arange(int(np.round(mn_elev-50, -2)), int(np.round(mx_elev+50, -2)), 100)

        total_area=pcr.cellvalue(pcr.maptotal(pcr.scalar(glacierShape)),1)[0]
        #print('total_area: '+str(total_area))

        if len(elevs)>1:
            E_norm=[(elevs[-1]-elevs[i])/(elevs[-1]-elevs[0]) for i in range(len(elevs))]

            #Following the values given by Huss et al., 2010
            if total_area>=20:
                a=-0.02
                g=6
                b=0.12
                c=0
            elif total_area<20 and total_area>=5:
                a=-0.05
                g=4
                b=0.19
                c=0.01
            elif total_area<5:
                a=-0.3
                g=2
                b=0.6
                c=0.09

            dH=[(E+a)**g+b*(E+a)+c for E in E_norm]

            #When dH is 0 and all other cells melt away, this can give an error. So needs to be a very small number.
            if dH[-1]==0:
                dH[-1]=1e-6

        else:
            dH=[1]


        total_diff=pcr.scalar(0.0)
        for k in range(len(fs_list)):#np.arange(0.01, 1, 0.01): 
            if k!=0:
                glacierShape=temp_volume>0
                volume=temp_volume

                assert abs((1-fs_list[k-1]/100)*M-pcr.cellvalue(pcr.maptotal(pcr.scalar(volume)),1)[0]) <= 1e-2

            glacierElev=pcr.ifthenelse(glacierShape, self.highResolutionDEM, np.nan)

            d=[]
            for j, elev in enumerate(elevs):
                area_dH=pcr.cellvalue(pcr.maptotal(pcr.scalar(pcr.pcrand(glacierElev>=float(elev), glacierElev<float(elev)+100))),1)[0]*dH[j]
                d+=[area_dH]

            Mnew=np.sum(d)



            if k!=0:
                #Every time remove the same amount of the mass from the glacier
                fs=(fs_list[1]*1e-2*M)/Mnew
            else:
                fs=0
            #print('fs: '+str(fs))

            fs_map=pcr.readmap("/hyclimm/gjanzing/data/glacierShape_minus"+str(int(fs_list[k])).zfill(2)+"M.map")            

            temp_volume=volume
            for j, elev in enumerate(elevs):
                temp_volume=pcr.ifthenelse(pcr.pcrand(glacierElev>=float(elev), glacierElev<(float(elev)+100)), temp_volume-fs*dH[j], temp_volume)

            #print('fs*dH: '+str(fs*np.sum(dH)))
            total_diff=pcr.ifthenelse(temp_volume<0.0, -1*temp_volume, 0.0)
            too_much_melt=pcr.cellvalue(pcr.maptotal(pcr.scalar(total_diff)),1)[0]
            temp_volume=pcr.ifthenelse(temp_volume<0.0, 0.0, temp_volume)
            #print('TMM: '+str(too_much_melt))


            #------FIXING WHEN THERE IS TOO MUCH MELT.....---------------------
            rounds=0
            while too_much_melt>=1e-2 or rounds==5:
                #print('Fixing: '+str(rounds))
                glacierShape=temp_volume>0

                glacierElev=pcr.ifthenelse(glacierShape, self.highResolutionDEM, np.nan)
                mn_elev = pcr.cellvalue(pcr.mapminimum(glacierElev),1)[0]
                mx_elev = pcr.cellvalue(pcr.mapmaximum(glacierElev),1)[0]
                elevs=np.arange(int(np.round(mn_elev-50, -2)), int(np.round(mx_elev+50, -2)), 100)

                d=[]
                for j, elev in enumerate(elevs):
                    area_dH=pcr.cellvalue(pcr.maptotal(pcr.scalar(pcr.pcrand(glacierElev>=float(elev), glacierElev<float(elev)+100))),1)[0]*dH[j]
                    d+=[area_dH]

                Mnew=np.sum(d)
                fs=too_much_melt/Mnew

                for j, elev in enumerate(elevs):
                    temp_volume=pcr.ifthenelse(pcr.pcrand(glacierElev>=float(elev), glacierElev<(float(elev)+100)), temp_volume-fs*dH[j], temp_volume)

                total_diff=pcr.ifthenelse(temp_volume<0.0, -1*temp_volume, 0.0)
                too_much_melt=pcr.cellvalue(pcr.maptotal(pcr.scalar(total_diff)),1)[0]
                temp_volume=pcr.ifthenelse(temp_volume<0.0, 0.0, temp_volume)
                rounds+=1

            if rounds==5:
                logger.info('ERROR----WOW NOT STABILISED YET!!!----ERROR-----ERROR----ERROR---')

            fs_map=fs_map+temp_volume
            pcr.report(fs_map,"/hyclimm/gjanzing/data/glacierShape_minus"+str(int(fs_list[k])).zfill(2)+"M.map")