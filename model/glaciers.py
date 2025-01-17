import netCDF4 as nc
import pcraster as pcr
import numpy as np

import logging
logger = logging.getLogger(__name__)

import virtualOS as vos
from ncConverter import *


def initializeGlacier(self, iniItems):
    """
    self: landSurface object
    iniItems: iniItems from landSurface module
    """
    
    logger.info("Initialising Glaciers...")
    self.glacierDir=iniItems.globalOptions['outputDir']

    #Print where glacier shapes are stored.
    print("----------GLACIER DIR:-----")
    print(self.glacierDir)

    #Load glacierIce
    self.glacierIce = vos.readPCRmapClone(\
                      iniItems.landSurfaceOptions['glaciers'],
                      self.cloneMap,self.tmpDir,self.inputDir)
    self.glacierIce=pcr.ifthenelse(self.glacierIce>0.0, self.glacierIce, 0.0)
    self.glacierIce_ini=self.glacierIce #Store initial glacier shapes for deltaH approach.
    self.glacierized=pcr.scalar(self.glacierIce>0.0) #Define glacierized terrain.
    self.glacierWater=pcr.scalar(0.0) #Initialize glacier water storage.

    #Load surface gradient.
    self.gradient = vos.readPCRmapClone(iniItems.routingOptions[str('gradient')],\
                             self.cloneMap,self.tmpDir,self.inputDir) 

    #Load glacierNumber (for deltaH approach) and make it scalar
    self.glacierNumber = vos.readPCRmapClone(\
                      iniItems.landSurfaceOptions['glacierNumber'],
                      self.cloneMap,self.tmpDir,self.inputDir)
    self.glacierNumber = pcr.scalar(self.glacierNumber)
    
    #Load LDDMap.
    self.lddMap = vos.readPCRmapClone(iniItems.routingOptions['lddMap'],\
                                            self.cloneMap,self.tmpDir,self.inputDir,True)
    self.lddMap = pcr.lddrepair(pcr.ldd(self.lddMap))
    self.lddMap = pcr.lddrepair(self.lddMap)

    
    #----Read Variables for updating glaciers:--------------------------------------------------------
    glacierParams      = ['degreeDayFactorGlacier','acclim','Kmin','Krange','Ag', 'freezingTGlacier'] #'AccEfficiency',

    for var in glacierParams:
        input = iniItems.landSurfaceOptions[str(var)]
        vars(self)[var] = vos.readPCRmapClone(input,self.cloneMap,
                                        self.tmpDir,self.inputDir)
        vars(self)[var] = pcr.spatial(pcr.scalar(vars(self)[var]))

    #Send glacier parameters to coverTypes.
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

        #In case we need to initialize deltaH (not necessary if same domain has been used before), do this.
        if self.ini_dH==True:
            self.highResolutionDEM = vos.readPCRmapClone(\
                                iniItems.meteoDownscalingOptions['highResolutionDEM'],
                                       self.cloneMap,self.tmpDir,self.inputDir)
            initializeDeltaH(self)
    
    if self.glacierType=='Immerzeel':

        #Parameters derived from Immerzeel et al., 2012
        self.tau_0=364000 #Equilibrium shear stress [Nm-2]
        self.R=5.8e10 #Material roughness coefficient [Nm-2s1/3]
        self.flowlim=5 #[m]; Flow limit. If ice is thinner, there is no flow
        
        #Calculate channel sizes and dimensions.
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
    
def updateStaticGlacier(self, meteo, currTimeStep):
    #Empirical constants
    #Glacier capacity
    glacierWaterHoldingCap=self.snowWaterHoldingCap
    #acclim=0 #minimum snow thickness needed for Accumulation [m]

    #Calculate melt rate of glacier by multiplying snow melt with a glacier factor; e.g. see Stahl et al., 2008; Seibert et al., 2018
    DDFGlacier=pcr.max(self.degreeDayFactor+self.degreeDayAmplitude*np.sin((currTimeStep.doy-81)*2*np.pi/366), 0.0)*self.degreeDayFactorGlacier

    #if currTimeStep.doy==30:
    #    a,b,c =vos.getMinMaxMean(pcr.spatial(DDFGlacier))
    #    msg = "Day:81: GlacierDDF initialized: Min %f Max %f Mean %f" %(a,b,c)
    #    logger.info(msg)
    #    print(msg)
    
    #Calculate Glacier Accumulation on before first of September.
    if currTimeStep.doy==244:
        #self.glacierAccumulation=pcr.ifthenelse(self.snowCoverSWE>self.acclim, self.AccEfficiency*self.glacierized*self.snowCoverSWE, 0)
        self.glacierAccumulation=pcr.ifthenelse(self.snowCoverSWE>self.acclim, self.snowCoverSWE-self.acclim, 0)
        self.glacierIce+=self.glacierAccumulation
        self.snowCoverSWE-=self.glacierAccumulation

        a,b,c =vos.getMinMaxMean(pcr.spatial(DDFGlacier))
        msg = "GlacierDDF initialized: Min %f Max %f Mean %f" %(a,b,c)
        logger.info(msg)
        print(msg)
    else:
        self.glacierAccumulation=pcr.scalar(0.0)

    #HBV START
    #For HBV
    #self.AccEfficiency=self.acclim
    #self.glacierAccumulation=pcr.ifthenelse(pcr.pcrand(self.glacierized>0.0, self.snowCoverSWE>0.0), self.AccEfficiency*self.snowCoverSWE, 0)
    #self.glacierIce+=self.glacierAccumulation
    #self.snowCoverSWE-=self.glacierAccumulation
    
        
    #Ice Melt
    #To Do: make melt different from snow; Add rain?
    #if self.seasonalMelt:
    #    self.iceMelt = \
    #        pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingT),(self.snowCoverSWE>0.0)), \
    #        0.0, \
    #       pcr.min(self.glacierIce, \
    #                self.glacierized*pcr.max(meteo.temperature - self.freezingT, 0.0)*(self.degreeDayFactor+self.degreeDayAmplitude*np.sin((currTimeStep.doy-81)*2*np.pi/366))))
    #else:

    #if self.freezingTGlacier==False:
    #    self.freezingTGlacier=self.freezingT        
    
    # self.iceMelt = \
    #     pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingTGlacier),(self.snowCoverSWE>0.0)), \
    #     0.0, \
    #    pcr.min(self.glacierIce, \
    #             self.glacierized*pcr.max(meteo.temperature - self.freezingTGlacier, 0.0) * (self.degreeDayFactorGlacier)))
    
    #self.iceMelt = \
    #    pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingTGlacier),(self.snowCoverSWE>0.0)), \
    #    0.0, \
    #   pcr.min(self.glacierIce, \
    #            self.glacierized*pcr.max(meteo.temperature - self.freezingTGlacier, 0.0) * (DDFGlacier)))

    #Calculate iceMelt using degree day approach.
    self.iceMelt = \
        pcr.ifthenelse(pcr.pcror((meteo.temperature <= self.freezingTGlacier),(self.snowCoverSWE>0.0)), \
        0.0, \
       pcr.min(self.glacierIce, \
                self.glacierized*pcr.max(meteo.temperature - self.freezingTGlacier, 0.0) * (DDFGlacier)))

    
    #Subtract iceMelt from glacier
    self.glacierIce-=self.iceMelt

    #For deltaH: keep track of yearly iceMelt and accumulation
    if self.glacierType=='delta_H':
        self.yearlyIceMelt+=self.iceMelt
        self.yearlyGlacierAcc+=self.glacierAccumulation
            
        
    #Update Glacier Water fluxes
    self.snow2GlacierWater= - pcr.min(0.0, self.deltaSnowCover) * self.glacierized #Convert to positive when added to glacier, remove refreezing
      
    self.rain2GlacierWater=pcr.ifthenelse(pcr.pcrand(self.glacierized>0.0, self.snowCoverSWE==0), self.liquidPrecip, 0.0)

    self.snowFreeWater = pcr.max(0.0, self.snowFreeWater \
                                 - self.deltaSnowCover - self.snow2GlacierWater \
                                 +  self.liquidPrecip -  self.rain2GlacierWater)
    
    self.capacity2GlacierWater=pcr.ifthenelse(self.glacierized>0.0, pcr.max(0., self.snowFreeWater - self.snowWaterHoldingCap * self.snowCoverSWE), 0.0)
   
    # update snowFreeWater (after capacity2GlacierWater) 
    self.snowFreeWater    = pcr.max(0., self.snowFreeWater - \
                                            self.capacity2GlacierWater)

    #Calculate glacierWater
    self.glacierWater=self.glacierWater+self.rain2GlacierWater+self.snow2GlacierWater+self.capacity2GlacierWater+self.iceMelt
    

    #Calculate Outflow, to be added to directRunoff.
    self.glacierOutflow=pcr.min(self.glacierWater, self.glacierWater*(self.Kmin+self.Krange*pcr.exp(-1*self.Ag*self.snowCoverSWE))*self.glacierized)
    self.glacierWater=pcr.max(0, self.glacierWater-self.glacierOutflow)


    #In case water there is too much water in the glacier (should not happen), create additional outflow
    self.netGlacierWaterExcess = pcr.max(0., self.glacierWater - \
                     glacierWaterHoldingCap * self.glacierIce)
    self.glacierWater=pcr.max(0, self.glacierWater-self.netGlacierWaterExcess)
    self.glacierOutflow=self.glacierOutflow+self.netGlacierWaterExcess

    

def glacierSlideImmerzeel(self, currTimeStep):
    #Calculate glacierSliding following Immerzeel et al., 2012.
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
    self.glacierIce=self.glacierIceVolume/self.cellArea
    
    self.glacierized=pcr.scalar(self.glacierIce>0.0)

    logger.info('Finished with GlacierSlideImmerzeel: Impressive....')

def updateDeltaH(self, currTimeStep):
    #Perform the deltaH calculation.
    logger.info('Updating delta H.......')

    #Initialize newGlacierIce and extraMelt variables.
    glacierIceNew=pcr.scalar(0.0)
    extraMelt=pcr.scalar(0.0)

    #Make the glacierNumber nominal.
    self.glacierNumber=pcr.nominal(self.glacierNumber)
    
    #Calculate Initial Volume per Glacier
    volume_ini=self.glacierIce_ini
    volume_ini=pcr.ifthenelse(volume_ini>0.0,volume_ini,0.0)
    M=pcr.areatotal(volume_ini, self.glacierNumber)

    #Calculate Volume at Present per Glacier
    volume=self.glacierIceYS
    volume=pcr.ifthenelse(volume>0.0,volume,0.0)
    Mpresent=pcr.areatotal(volume, self.glacierNumber)

    #Store previous glacierNumber
    oldGlacierNumber=pcr.ifthenelse(volume>0,self.glacierNumber,0)


    #Make variables spatial, if they are not already.
    self.yearlyGlacierAcc=pcr.spatial(self.yearlyGlacierAcc)
    self.yearlyIceMelt=pcr.spatial(self.yearlyIceMelt)

    #Calculate Mass Balance over the Year per Glacier
    self.glacierMB=self.yearlyGlacierAcc-self.yearlyIceMelt    
    deltaM=pcr.areatotal(self.glacierMB,oldGlacierNumber)

    #Calculate the mass fraction of the new glacier.
    frac_min=(deltaM+Mpresent)/M

    #Print in between information.
    a,b,c =vos.getMinMaxMean(pcr.spatial(self.glacierMB))
    msg = "GlacierMB initialized: Min %f Max %f Mean %f" %(a,b,c)
    logger.info(msg)
    print(msg)

    a,b,c =vos.getMinMaxMean(pcr.spatial(self.yearlyIceMelt))
    msg = "Glacier yearlyIceMelt initialized: Min %f Max %f Mean %f" %(a,b,c)
    logger.info(msg)
    print(msg)

    a,b,c =vos.getMinMaxMean(pcr.spatial(self.yearlyGlacierAcc))
    msg = "self.yearlyGlacierAcc initialized: Min %f Max %f Mean %f" %(a,b,c)
    logger.info(msg)
    print(msg)

    a,b,c =vos.getMinMaxMean(pcr.spatial(frac_min))
    msg = "Frac_min initialized: Min %f Max %f Mean %f" %(a,b,c)
    logger.info(msg)
    print(msg)

    a,b,c =vos.getMinMaxMean(pcr.spatial(deltaM))
    msg = "deltaM initialized: Min %f Max %f Mean %f" %(a,b,c)
    logger.info(msg)
    print(msg)

    a,b,c =vos.getMinMaxMean(pcr.spatial(Mpresent))
    msg = "Mpresent initialized: Min %f Max %f Mean %f" %(a,b,c)
    logger.info(msg)
    print(msg)

    a,b,c =vos.getMinMaxMean(pcr.spatial(M))
    msg = "M initialized: Min %f Max %f Mean %f" %(a,b,c)
    logger.info(msg)
    print(msg)

    
    #Calculate fraction of mass loss
    frac=1-frac_min
    frac=pcr.ifthenelse(frac>1, 1, frac)
    frac=pcr.ifthenelse(frac<0, 0, frac)
    frac=pcr.ordinal(pcr.roundoff(frac*100))
    
    
    #Load all the new glacier shapes.
    #This is one loop that I can't really avoid :(
    for i in range(101):
        shapes=pcr.readmap(self.glacierDir+"glacierShape_minus"+str(int(i)).zfill(2)+"M.map")
        newShape=pcr.ifthenelse(pcr.ordinal(frac)==pcr.ordinal(i), shapes, 0.0)
        newShape=pcr.ifthenelse(newShape>0, newShape, 0.0)
        newShape=pcr.cover(newShape, 0.0)
        glacierIceNew=glacierIceNew+newShape #Update the new glacier shapes

    #We don't want to have glacier growing when they are actually loosing mass due to round-off errors.
    #Therefore, in case a glacier does not move to another percentage point, we just take the old shape
    old_frac=1-(Mpresent)/M
    old_frac=pcr.ifthenelse(old_frac>1, 1, old_frac)
    old_frac=pcr.ifthenelse(old_frac<0, 0, old_frac)
    old_frac=pcr.ordinal(pcr.roundoff(old_frac*100))            
    glacierIceNew=pcr.ifthenelse(old_frac==frac, self.glacierIceYS, glacierIceNew) 

    #Now: Fixing the water balance! Redistribute any potential mismatches in the water balance.
    Mnew=pcr.areatotal(glacierIceNew, self.glacierNumber)
    glacierGap=(Mpresent-Mnew)+deltaM #Gap in mass balance

    #If difference between previous shape and new shape is larger than delta M (glacierDiff is positive),
    #this means that there is more melt than is observed.
    #Extra melt has to be added to the glacier.
    #And Vice Versa
    
    #Only divide too much or too little melt over the glacier.
    #NOTE: also when a glacier grows more than the initial shape, it will be divided over the glacier as glacier ice!!

    #We only divide stuff based on the current thickness, in order to avoid melting too much in cells with hardly any ice.
    #pcr.setglobaloption("unitcell") #--> this might have consequences?
    ratio=glacierIceNew/Mnew #We need to divide it only over the area where there is now glacier.
    
    glacierGap=pcr.ifthenelse(glacierIceNew>0, glacierGap, 0) #We need to divide it only over the area where there is now glacier.
    glacierDiff=pcr.ifthenelse(pcr.pcrand(pcr.ordinal(frac)!=100, pcr.abs(glacierGap)>1e-3), pcr.scalar(glacierGap*ratio), 0.0) #Division based on ice thickness.
    
    glacierDiff=pcr.cover(glacierDiff, 0.0) #Make sure no NaNs exist (due to ratio?)
    glacierGap=pcr.cover(glacierGap, 0.0) #Make sure no NaNs exist (due to ratio?)
    frac=pcr.cover(frac, 0) #Make sure no NaNs exist (due to ratio?)

    #Glaciers can not grow beyond their original extend.
    #Transport additional accumulation (larger than initial state) downstream.
    positiveGap=pcr.ifthenelse(pcr.pcror((frac<0), pcr.pcrand((frac==0),(glacierGap>0))), glacierDiff, 0.0)
    glacierDiff=pcr.ifthenelse(pcr.pcror((frac<0), pcr.pcrand((frac==0),(glacierGap>0))), 0, glacierDiff)
    
    glacierDiff=pcr.ifthenelse(pcr.pcror((frac<0), pcr.pcrand((frac==0),(glacierGap>0))), 0, glacierDiff)
    transportFraction=pcr.ifthenelse(self.glacierIce>0, pcr.scalar(1), pcr.scalar(0)) #All ice is transported downwards
    added = pcr.accufractionstate(self.lddMap, positiveGap*self.cellArea, transportFraction)
    added=added/self.cellArea
    glacierDiff=glacierDiff+added
    glacierIceNew=glacierIceNew+glacierDiff #Update new glacierShapes
     
    #Dealing with too much melt: Should not happen: 
    #Melt is restricted to the glacier volume.
    extraMelt=pcr.min(glacierIceNew, 0.0)
    glacierIceNew=pcr.max(glacierIceNew, 0.0)
        
    #Initial quick water balance check:
    #assert pcr.cellvalue(pcr.mapminimum(pcr.scalar(extraMelt)),1)[0] >= -0.01
    extraMelt=pcr.scalar(extraMelt)
    em_error = pcr.cellvalue(pcr.mapminimum(extraMelt),1)[0]
    a,b,c =vos.getMinMaxMean(extraMelt)
    msg = "Extra melt: Min %f Max %f Mean %f" %(a,b,c)
    print(msg)
    if em_error<-1e-4:
        logger.info('ERROR: EXTRA MELT: '+str(em_error))

    #Update the new GlacierIce!
    self.glacierIce=glacierIceNew 
    #FIX: cover NaN values
    self.glacierIce=pcr.cover(self.glacierIce, 0.0)
    self.glacierized=pcr.scalar(self.glacierIce>0.0)

    #Re-initialize the storages over time.
    self.glacierIceYS=self.glacierIce
    self.yearlyIceMelt=pcr.scalar(0.0)
    self.yearlyGlacierAcc=pcr.scalar(0.0)
    logger.info('Finished with delta H!!.......')
    


# def updateDeltaH_old(self, currTimeStep):
#     logger.info('Updating delta H.......')
#     mn = pcr.cellvalue(pcr.mapminimum(self.glacierNumber),1)[0]
#     mx = pcr.cellvalue(pcr.mapmaximum(self.glacierNumber),1)[0]

#     glacierIcenew=pcr.scalar(0.0)
#     extraMelt=pcr.scalar(0.0)
    
#     self.glacierMB=self.yearlyGlacierAcc-self.yearlyIceMelt

#     for i in range(int(mn), int(mx)+1):
#     #for i in range(int(mn), int(mn)+30):
#         glacierShape=self.glacierNumber==i
#         volume_ini=pcr.ifthenelse(glacierShape, self.glacierIce_ini, 0.0)
#         volume_ini=pcr.ifthenelse(volume_ini>0.0,volume_ini,0.0)
#         glacierShape=volume_ini>0
#         valid=pcr.cellvalue(pcr.maptotal(pcr.scalar(volume_ini>0)),1)[0]

#         if valid==0:
#             continue

#         M=pcr.cellvalue(pcr.maptotal(volume_ini),1)[0]

#         volume=pcr.ifthenelse(glacierShape, self.glacierIceYS, 0.0)
#         volume=pcr.ifthenelse(volume>0.0,volume,0.0)
#         Mpresent=pcr.cellvalue(pcr.maptotal(volume),1)[0]

#         indivMB=pcr.ifthenelse(glacierShape, self.glacierMB, 0.0)
#         deltaM=pcr.cellvalue(pcr.maptotal(indivMB),1)[0]

#         frac_min=np.round((deltaM+Mpresent)/M,2)
#         frac=1-frac_min
#         if frac>1:
#             frac=1

#         #We don't want to have glacier growing when they are actually loosing mass due to round-off errors.
#         #Therefore, in case a glacier does not move to another percentage point, we just take the old shape
#         old_frac=1-np.round((Mpresent)/M,2)
#         if int(frac*100)==int(old_frac*100):
#             indivNew=self.glacierIceYS
#         elif frac>0:
#             indivNew=pcr.readmap("/hyclimm/gjanzing/data/glacierShape_minus"+str(int(frac*100)).zfill(2)+"M.map")
#         elif frac<0:
#             #print('Sir, we have caught a grower!!')
#             indivNew=pcr.readmap("/hyclimm/gjanzing/data/glacierShape_minus00M.map")
        
#         #FIX: cover NaN values
#         indivNew=pcr.cover(indivNew, 0.0)
        
#         newShape=pcr.ifthenelse(glacierShape, indivNew, 0.0)
#         glacierIcenew=glacierIcenew+newShape

#         Mnew=pcr.cellvalue(pcr.maptotal(newShape),1)[0]

#         #To Do: what to do with extra volume??
#         #--> add to of subtract from glacier Ice
#         # else groundwater?

#         #If difference between previous shape and new shape is larger than delta M (glacierDiff is positive),
#         #this means that there is more melt than is observed.
#         #Extra melt has to be added to the glacier.
#         #And Vice Versa
#         glacierGap=(Mpresent-Mnew)+deltaM

#         #print('Mnew: '+str(Mnew))
#         #print('actual delta: '+str(Mnew-Mpresent))

#         #Only divide too much or too little melt over the glacier.
#         #NOTE: also when a glacier grows more than the initial shape, it will be divided over the glacier as glacier ice!!

#         if (frac!=1) and (abs(glacierGap)>1e-3):
#             #We only divide stuff over places where the thickness is larger than 1, in order to prevent cells where the glacierDiff is larger than the ice thickness.
#             area=pcr.cellvalue(pcr.maptotal(pcr.scalar(newShape>1)),1)[0]

#             #If all cells are smaller than 1m, than area=1.
#             if area==0:
#                 area=1

#             #glacierDiff=pcr.ifthenelse(glacierShape, pcr.scalar(glacierGap/area), 0.0)
#             glacierDiff=pcr.ifthenelse(newShape>1, pcr.scalar(glacierGap/area), 0.0)
            
#             #Glaciers can not grow beyond their original extend. In case this happens, it is added to snow. 
#             #This can then lead to snow towers, which could be removed by snow redistribution.
#             #if (frac<0) or ((frac==0) and (glacierGap>0)):
#             #    self.snowCoverSWE+=glacierDiff
                
#             #    #This will be overwritten anyway, right?
#             #    self.glacierAccuCorrection+=glacierDiff
#             #else:
#             glacierIcenew=glacierIcenew+glacierDiff
        
        
#         #Should not happen: melt is restricted to the glacier volume, right?
#         #Should this be though??
#         extraMelt+=pcr.min(glacierIcenew, 0.0)
#         glacierIcenew=pcr.max(glacierIcenew, 0.0)
        
#         #Implement a good water balance check:..... but how?
#         #assert pcr.cellvalue(pcr.mapminimum(pcr.scalar(extraMelt)),1)[0] >= -0.01
    
#     em_error = pcr.cellvalue(pcr.mapminimum(extraMelt),1)[0]
    
#     a,b,c =vos.getMinMaxMean(extraMelt)
#     msg = "Extra melt: Min %f Max %f Mean %f" %(a,b,c)
#     print(msg)

#     if em_error<-1e-4:
#         logger.info('ERROR: EXTRA MELT: '+str(em_error))

#     #Update the new GlacierIce!
#     self.glacierIce=glacierIcenew 
#     #FIX: cover NaN values
#     self.glacierIce=pcr.cover(self.glacierIce, 0.0)
#     self.glacierized=pcr.scalar(self.glacierIce>0.0)

#     #Re-initialize the storages over time.
#     self.glacierIceYS=self.glacierIce
#     self.yearlyIceMelt=pcr.scalar(0.0)
#     self.yearlyGlacierAcc=pcr.scalar(0.0)
#     logger.info('Finished with delta H!!.......')
def initializeDeltaH(self):
    """
    Initialize deltaH approach.
    """
    
    logger.info('Initializing Delta H parameterization....')
    print('Initializing Delta H parameterization....')
    
    mn = pcr.cellvalue(pcr.mapminimum(self.glacierNumber),1)[0]
    mn=np.max([mn, 1])
    mx = pcr.cellvalue(pcr.mapmaximum(self.glacierNumber),1)[0]
    self.glacierNumber=pcr.nominal(self.glacierNumber)

    #Create glaciers shapes with a certain mass loss and stores these.
    fs_list=np.arange(0, 101, 1)
    for frac in fs_list:
        print(frac)
        fs_map=pcr.scalar(0.0)
        logger.info("Creating: " +self.glacierDir+"glacierShape_minus"+str(int(frac)).zfill(2)+"M.map")
        pcr.report(fs_map, self.glacierDir+"glacierShape_minus"+str(int(frac)).zfill(2)+"M.map")

    #Calculate Initial Volume per Glacier
    volume=self.glacierIce
    volume=pcr.ifthenelse(volume>0.0,volume,0.0)
    M=pcr.areatotal(volume, self.glacierNumber)
    
    #Calculate Elevations
    glacierElev=pcr.ifthenelse(self.glacierIce>0, self.highResolutionDEM, 0)
    mn_elev = pcr.areaminimum(glacierElev, self.glacierNumber)
    mx_elev = pcr.areamaximum(glacierElev, self.glacierNumber)

    #elev_steps=pcr.roundoff((mx_elev+50)/100)-pcr.roundoff((mn_elev-50)/100)
    elev={}
    steps=20 #Use 20 elevation bands per glacier.
    deltaElev=(pcr.roundoff((mx_elev+50)/100)-pcr.roundoff((mn_elev-50)/100))/steps*100 #Round to 100 and then divide into 20 steps. #Shouldn't divide into hundred steps, but into steps of 100....
    for i in range(steps):
        elev[i]=pcr.roundoff((mn_elev-50)/100)*100+i*deltaElev

    #Calculate total surface area per glacier(in number of cells)
    pcr.setglobaloption("unitcell")
    total_area=pcr.areaarea(self.glacierNumber)
    pcr.setglobaloption("unittrue")

    #Normalize elevation zones.
    E_norm={}
    for i in range(steps):
        E_norm[i]=(elev[steps-1]-elev[i])/(elev[steps-1]-elev[0])

    #Calculate shape parameters of glaciers
    a=pcr.scalar(0.0)
    g=pcr.scalar(0.0)
    b=pcr.scalar(0.0)
    c=pcr.scalar(0.0)

    a=pcr.ifthenelse(total_area>=20,-0.02,a)
    a=pcr.ifthenelse(pcr.pcrand(total_area<20,total_area>=5),-0.05,a)
    a=pcr.ifthenelse(total_area<5,-0.3,a)

    g=pcr.ifthenelse(total_area>=20,6,g)
    g=pcr.ifthenelse(pcr.pcrand(total_area<20,total_area>=5),4,g)
    g=pcr.ifthenelse(total_area<5,2,g)

    b=pcr.ifthenelse(total_area>=20,0.12,b)
    b=pcr.ifthenelse(pcr.pcrand(total_area<20,total_area>=5),0.19,b)
    b=pcr.ifthenelse(total_area<5,0.6,b)

    c=pcr.ifthenelse(total_area>=20,0,c)
    c=pcr.ifthenelse(pcr.pcrand(total_area<20,total_area>=5),0.01,c)
    c=pcr.ifthenelse(total_area<5,0.09,c)

    #Calculate the dH function.
    dH={}
    for i in range(steps):
        dH[i]=(E_norm[i]+a)**g+b*(E_norm[i]+a)+c
        #When dH is 0 and all other cells melt away, this can give an error. So needs to be a very small number.
        dH[i]=pcr.ifthenelse(dH[i]==0, 1e-6, dH[i])

    total_diff=pcr.scalar(0.0)

    #Now, create maps for each percentage massloss.
    for k in range(len(fs_list)):
        print("Step:"+str(k))
        if k!=0:
            volume=temp_volume #Volume refers to the previous time steps.
            
            #Perform some checks for debugging.
            prevVol=pcr.areatotal(volume, self.glacierNumber)
            x,y,z =vos.getMinMaxMean(prevVol)
            msg = "Volume: Min %f Max %f Mean %f" %(x,y,z)
            print(msg)

            #Calculate difference between new mass and what it should be?
            iniMassDiff=(1-fs_list[k-1]/100)*M-pcr.areatotal(volume, self.glacierNumber)

            #Perform some checks for debugging: is there an error in the initial mass of the glaciers?
            x,y,z =vos.getMinMaxMean(iniMassDiff)
            msg = "Ini Mass Difference: Min %f Max %f Mean %f" %(x,y,z)
            print(msg)
            #assert iniMassDiff <= 1e-1

        #Calculate total normalized mass changed (including a correction for the surface area per elevation zone)
        glacierElev=pcr.ifthenelse(volume>0, self.highResolutionDEM, 0) #glacier elevation
        area_dH={}
        Mnew=pcr.scalar(0.0)
        for i in range(steps):
            if i<(steps-1):
                elev_range=pcr.pcrand(glacierElev>=elev[i], glacierElev<elev[i+1])
            else:
                elev_range=pcr.pcrand(glacierElev>=elev[i], glacierElev<elev[i]+100)
            area_dH[i]=pcr.areatotal(pcr.scalar(elev_range), self.glacierNumber)*dH[i] #Correct dH for the surface area of this elevation zone.
            Mnew=Mnew+area_dH[i] #Summed the normalized mass changes

        if k!=0:
            #Every time remove the same amount of the mass from the glacier: e.g. a percentage from the original volume.
            fs=(fs_list[1]*1e-2*M)/Mnew #Calculate the scaling factor
        else:
            fs=0

        temp_volume=volume
        
        #Can't merge this loop because we need to know the new volume to calculate fs...; or can I?
        for i in range(steps):
            if i<(steps-1):
                elev_range=pcr.pcrand(glacierElev>=elev[i], glacierElev<elev[i+1])
            else:
                elev_range=pcr.pcrand(glacierElev>=elev[i], glacierElev<elev[i]+100)
            temp_volume=pcr.ifthenelse(elev_range, temp_volume-fs*dH[i], temp_volume) #Calculate the new glacier volume and distribution by applying the scaling factor to the deltaH function.
            
        #Adjust for small volume errors.
        total_diff=pcr.ifthenelse(temp_volume<0.0, -1*temp_volume, 0.0)
        too_much_melt=pcr.areatotal(total_diff, self.glacierNumber)
        temp_volume=pcr.ifthenelse(temp_volume<0.0, 0.0, temp_volume)

        ratio=temp_volume/pcr.areatotal(temp_volume, self.glacierNumber)
        temp_volume=temp_volume-too_much_melt*ratio

        #Store the final map.
        fs_map=temp_volume
        fs_map=pcr.cover(fs_map, 0.0)
        pcr.report(fs_map, self.glacierDir+"glacierShape_minus"+str(int(k)).zfill(2)+"M.map")
            
# def initializeDeltaH(self):
#     logger.info('Initializing Delta H parameterization....')

#     mn = pcr.cellvalue(pcr.mapminimum(self.glacierNumber),1)[0]
#     mx = pcr.cellvalue(pcr.mapmaximum(self.glacierNumber),1)[0]


#     #fs_list=np.arange(0.01, 1, 0.01)
#     fs_list=np.arange(0, 101, 1)
#     #fs_list=np.arange(0, 1+1e-5, 0.1)
#     for frac in fs_list:
#         print(frac)
#         fs_map=pcr.scalar(0.0)
#         logger.info("Creating: /hyclimm/gjanzing/data/glacierShape_minus"+str(int(frac)).zfill(2)+"M.map")
#         pcr.report(fs_map,"/hyclimm/gjanzing/data/glacierShape_minus"+str(int(frac)).zfill(2)+"M.map")

#     #mn=16823
#     #print(mx)
#     for i in range(int(mn), int(mx)+1):

#         glacierShape=self.glacierNumber==i
#         volume=pcr.ifthenelse(glacierShape, self.glacierIce, 0.0)
#         volume=pcr.ifthenelse(volume>0.0,volume,0.0)
#         glacierShape=volume>0
#         valid=pcr.cellvalue(pcr.maptotal(pcr.scalar(volume>0)),1)[0]

#         if valid==0:
#             #When there are no cells, continue
#             continue

#         M=pcr.cellvalue(pcr.maptotal(volume),1)[0]

#         glacierElev=pcr.ifthenelse(glacierShape, self.highResolutionDEM, np.nan)
#         mn_elev = pcr.cellvalue(pcr.mapminimum(glacierElev),1)[0]
#         mx_elev = pcr.cellvalue(pcr.mapmaximum(glacierElev),1)[0]

#         elevs=np.arange(int(np.round(mn_elev-50, -2)), int(np.round(mx_elev+50, -2)), 100)

#         total_area=pcr.cellvalue(pcr.maptotal(pcr.scalar(glacierShape)),1)[0]
#         #print('total_area: '+str(total_area))

#         if len(elevs)>1:
#             E_norm=[(elevs[-1]-elevs[i])/(elevs[-1]-elevs[0]) for i in range(len(elevs))]

#             #Following the values given by Huss et al., 2010
#             if total_area>=20:
#                 a=-0.02
#                 g=6
#                 b=0.12
#                 c=0
#             elif total_area<20 and total_area>=5:
#                 a=-0.05
#                 g=4
#                 b=0.19
#                 c=0.01
#             elif total_area<5:
#                 a=-0.3
#                 g=2
#                 b=0.6
#                 c=0.09

#             dH=[(E+a)**g+b*(E+a)+c for E in E_norm]

#             #When dH is 0 and all other cells melt away, this can give an error. So needs to be a very small number.
#             if dH[-1]==0:
#                 dH[-1]=1e-6

#         else:
#             dH=[1]


#         total_diff=pcr.scalar(0.0)
#         for k in range(len(fs_list)):#np.arange(0.01, 1, 0.01): 
#             if k!=0:
#                 glacierShape=temp_volume>0
#                 volume=temp_volume

#                 assert abs((1-fs_list[k-1]/100)*M-pcr.cellvalue(pcr.maptotal(pcr.scalar(volume)),1)[0]) <= 1e-2

#             glacierElev=pcr.ifthenelse(glacierShape, self.highResolutionDEM, np.nan)

#             d=[]
#             for j, elev in enumerate(elevs):
#                 area_dH=pcr.cellvalue(pcr.maptotal(pcr.scalar(pcr.pcrand(glacierElev>=float(elev), glacierElev<float(elev)+100))),1)[0]*dH[j]
#                 d+=[area_dH]

#             Mnew=np.sum(d)



#             if k!=0:
#                 #Every time remove the same amount of the mass from the glacier
#                 fs=(fs_list[1]*1e-2*M)/Mnew
#             else:
#                 fs=0
#             #print('fs: '+str(fs))

#             fs_map=pcr.readmap("/hyclimm/gjanzing/data/glacierShape_minus"+str(int(fs_list[k])).zfill(2)+"M.map")            

#             temp_volume=volume
#             for j, elev in enumerate(elevs):
#                 temp_volume=pcr.ifthenelse(pcr.pcrand(glacierElev>=float(elev), glacierElev<(float(elev)+100)), temp_volume-fs*dH[j], temp_volume)

#             #print('fs*dH: '+str(fs*np.sum(dH)))
#             total_diff=pcr.ifthenelse(temp_volume<0.0, -1*temp_volume, 0.0)
#             too_much_melt=pcr.cellvalue(pcr.maptotal(pcr.scalar(total_diff)),1)[0]
#             temp_volume=pcr.ifthenelse(temp_volume<0.0, 0.0, temp_volume)
#             #print('TMM: '+str(too_much_melt))


#             #------FIXING WHEN THERE IS TOO MUCH MELT.....---------------------
#             rounds=0
#             while too_much_melt>=1e-2 or rounds==5:
#                 #print('Fixing: '+str(rounds))
#                 glacierShape=temp_volume>0

#                 glacierElev=pcr.ifthenelse(glacierShape, self.highResolutionDEM, np.nan)
#                 mn_elev = pcr.cellvalue(pcr.mapminimum(glacierElev),1)[0]
#                 mx_elev = pcr.cellvalue(pcr.mapmaximum(glacierElev),1)[0]
#                 elevs=np.arange(int(np.round(mn_elev-50, -2)), int(np.round(mx_elev+50, -2)), 100)

#                 d=[]
#                 for j, elev in enumerate(elevs):
#                     area_dH=pcr.cellvalue(pcr.maptotal(pcr.scalar(pcr.pcrand(glacierElev>=float(elev), glacierElev<float(elev)+100))),1)[0]*dH[j]
#                     d+=[area_dH]

#                 Mnew=np.sum(d)
#                 fs=too_much_melt/Mnew

#                 for j, elev in enumerate(elevs):
#                     temp_volume=pcr.ifthenelse(pcr.pcrand(glacierElev>=float(elev), glacierElev<(float(elev)+100)), temp_volume-fs*dH[j], temp_volume)

#                 total_diff=pcr.ifthenelse(temp_volume<0.0, -1*temp_volume, 0.0)
#                 too_much_melt=pcr.cellvalue(pcr.maptotal(pcr.scalar(total_diff)),1)[0]
#                 temp_volume=pcr.ifthenelse(temp_volume<0.0, 0.0, temp_volume)
#                 rounds+=1

#             if rounds==5:
#                 logger.info('ERROR----WOW NOT STABILISED YET!!!----ERROR-----ERROR----ERROR---')

#             fs_map=fs_map+temp_volume
#             pcr.report(fs_map,"/hyclimm/gjanzing/data/glacierShape_minus"+str(int(fs_list[k])).zfill(2)+"M.map")