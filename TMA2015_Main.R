# Thermal Model for Anaerobic Digesters
# R_TMAD v.0.1
# Written in R version 3.4.0
# Last update: 10/05/2017
# Simon V. Pedersen (May 2017)


####################################################################################################################
#                                                       CONSTANTS
####################################################################################################################

# Physical, Chemical and atmospheric variables considered constant in the model
sigma <- 5.67037*10^(-8)               # Stefan-Boltzmann Constant                               (W m^-2 K^-4) , CODATA 2014
P0 <- 101325                           # Standard Atmosphere. Pressure at sea-level              (Pa)          , CODATA 2014
Sp0 <- 1360                            # Extraterretrial flux density normal to the solar beam   (W m^-2)      , Campbell 1977 

# Constant heat transfer coefficients
h.cov_air <- 3.55                      # Convective HTC Cover to Ambient Air                     (W m^-2 K^-1)
h.cov_gas <- 2.15                      # Convective HTC cover to biogas                          (W m^-2 K^-1)
h.gas_wall <- 2.70                     # Convective HTC Biogas to wall (headspace)               (W m^-2 K^-1)
h.gas_sub <- 2.20                      # Convective HTC Biogas to substrate                      (W m^-2 K^-1)
h.sub_wall <- 177.25                   # Convective HTC Substrate to wall                        (W m^-2 K^-1)
h.sub_floor <- 244.45                  # Convective HTC Substrate to floor                       (W m^-2 K^-1)


####################################################################################################################
#                                               MODEL COMPUTATION ENGINE
####################################################################################################################


##################################################################
# Function: OverallHeatTransferCoef
# Purpose: Calculates the overall heat transfer coefficient 
# Prog. Id.: SVIP
# Date: 15/05/2017
# Returns: The overll heat transfer coefficient (W m^-2 K^-1)
##################################################################

TMAMModel <- function(
  timeZone,                                 # Location time-zone from GMT                                                  (hr)
  longitude,                                # Location longitude                                                           (DEGREES)
  latitude,                                 # Location latitude                                                            (DEGREES)                                  
  masl,                                     # Digester location above sea-level                                            (m)
  aCoef,                                    # Coefficient of Tranismissivity. Constant. Average annual condition           (Dimensionless)
  dosStop,                                  # Number of days to conduct the simulation                                     (Integer)
  startDOY,                                 # The Day of Year at which the simulation starts                               (Integer)
  timeStep,                                 # The time-step for the Euler Forward method. Divisible by 1 hour (60*60 s)    (s)
  tempSubstrate0,                           # Initial temperature of the substrate                                         (K)
  annualAverageTemperature,                 # Annual ambient air temperature average                                       (K)
  annualAverageTemperatureAmplitude,        # Annual ambient air temperature amplitude                                     (K)
  coldestDOY,                               # Coldest Day of Year                                                          (Integer)
  annAvgSurfTemperature,                    # Annual average soil surface temperature                                      (K)
  annAvgSurfTempAmplitude,                  # Annual soil surface temperature amplitude                                    (K)
  wallThickness,                            # Digester wall thickness                                                      (m)
  radiusCover,                              # Combined radius of the cover (of the two covers)                             (m)
  radiusSubstrateInterface,                 # Radius of the substrate-biogas interface                                     (m)
  volumeSubstrate,                          # Constant volume of the substrate in the digester                             (m^3)
  digesterHeight,                           # The total height of the digester                                             (m)
  substrateDepth,                           # The depth of the substrate from bottom to surface                            (m)
  surfaceAreaWallBiogas,                    # Surface area of the wall-biogas interface                                    (m^2)
  surfaceAreaWallSubstrate,                 # Surface area of the wall-substrate interface                                 (m^2)
  surfaceAreaWallFloor,                     # Surface area of the floor-substrate interface                                (m^2)
  subtrateLoadingRate,                      # Constant substrate loading rate                                              (kg s^-1)
  kDigester,                                # Thermal conductivity of the digester material                                (W m^-1 k^-1)
  rhoSubstrate,                             # Density of the substrate. Assumed equivalent to water.                       (kg m^-3)
  cpSubstrate,                              # Specific heat of the substrate (slurry)                                      (J kg^-1 K^-1)
  alphaSoil,                                # Thermal diffusivity of the soil                                              (m^2 s^-1)
  emisCover,                                # Emissivity of the cover                                                      (DIMENSIONLESS)
  emisSubstrate,                            # Emissivity of the substrate                                                  (DIMENSIONLESS)
  abCover,                                  # Absorptivity of the cover                                                    (DIMENSIONLESS)
  heatingRate                               # External heating rate                                                        (W)
){

  
  ############################################################################
  # INITIALIZATION OF MODELS, PARAMETERS, ETC.
  ############################################################################
  # Heat transfer parameters
  frameResistanceBiogas <- data.frame(deltaL=c(wallThickness,1,1),kh=c(kDigester,h.gas_wall,h.gas_sub))
  frameResistanceSubstrate <- data.frame(deltaL=c(wallThickness,1),kh=c(kDigester,h.sub_wall))
  frameResistanceFloor <- data.frame(deltaL=c(wallThickness,1),kh=c(kDigester,h.sub_floor))
  frameResistanceCover <- data.frame(deltaL=c(1,wallThickness,1,1),kh=c(h.cov_air,kDigester,h.cov_gas,h.gas_sub))
  
  # Soil descretization frames
  descrecUnit <- 10
  seqDescretizationBiogas <- seq(0,(digesterHeight-substrateDepth),(digesterHeight-substrateDepth)/descrecUnit)
  seqDescretizationSubstrate <- seq((digesterHeight-substrateDepth),digesterHeight,(digesterHeight-(digesterHeight-substrateDepth))/descrecUnit)
  
  # Calculate Overall Heat Transfer Coefficients
  uOverallBiogas <- OverallHeatTransferCoef(frameResistanceBiogas)
  uOverallSubstrate <- OverallHeatTransferCoef(frameResistanceSubstrate)
  uOverallFloor <- OverallHeatTransferCoef(frameResistanceFloor)
  uOverallCover <- OverallHeatTransferCoef(frameResistanceCover)
  
  areaCover <- (radiusCover^2)*pi                             # Area of the digester cover                                                   (m^2)
  surfaceAreaSubstrate <- (radiusSubstrateInterface^2)*pi     # Surface area of the biogas-substrate interface                               (m^2)
  heightHeadspace <- (digesterHeight - substrateDepth)        # The total height of the digester headspace                                   (m)
  tempSubstrate <- tempSubstrate0
  substrateLoadingRate <- 0
  
  #message("areaCover:" , areaCover)
  #message("surfaceAreaSubstrate:" , surfaceAreaSubstrate)
  #message("heightHeadspace:" , heightHeadspace)
  #message("tempSubstrate:" , tempSubstrate)
  #message("substrateLoadingRate:" , substrateLoadingRate)
  
  #message("uOverallBiogas: ", uOverallBiogas)
  #message("uOverallSubstrate: ", uOverallSubstrate)
  #message("uOverallFloor: ", uOverallFloor)
  #message("uOverallCover: ", uOverallCover)

  
  #message((((60*60)/timeStep)*24*dosStop))
  # Initialize result matrix
  mat.res <- matrix(0,nrow=(((60*60)/timeStep)*24*dosStop),ncol=(5+7+3),dimnames = list(NULL,c("DOS","Hour","DOYTime","AirTemp","SubstrateTemp","qIRRa4","qRADa4","qCON42gr","qCON4gr","qCONfloorgr","qCONa4","qADVload","SoilTempSub","SoilTempGas","SoilTempFloor")))
  
  # Loop control parameters
  i <- 1
  
  # Geometric setup
  viewFactor.cov_sub <- CylendricalViewFactor(radiusCover,radiusSubstrateInterface, heightHeadspace)
  #message("ViewFactor:" , viewFactor.cov_sub)
  
  # Soil setup
  soilDampingDepthYr <- DampingDepth(((2*pi)/(365*24*60*60)),alphaSoil)
  #message("DampingDepthYr:" , soilDampingDepthYr)

  ############################################################################
  # SIMULATION LOOP
  ############################################################################
  for(dos in 1:dosStop){
    
    hr <- 0
  
    while(hr < 23.999){
      
      # Calculate the local time as a decimal DOY, at the current timestep
      localDOYTime <- startDOY + (dos-1) + (hr/24)
      totalDOY <- startDOY + (dos-1)
      currentDOY <- 0
      if(totalDOY > 365){
        currentDOY <- 1
      }else{
        currentDOY <- totalDOY
      }
      
      # Calculate the ambient air temperature at the current time-step, and other temperatures
      tempAir <- AmbientAirTemperature(annualAverageTemperature,annualAverageTemperatureAmplitude,coldestDOY,localDOYTime)
      effSkyTemperature <- SkyTemperature(tempAir)
      tempLoading <- tempAir + 1
      
      # Calculate mean soil temperature by the descretization
      meanSoilTempSub   <- mean(SoilTemperature(seqDescretizationSubstrate,localDOYTime,annAvgSurfTemperature,annAvgSurfTempAmplitude,(2*pi/(365)),soilDampingDepthYr,coldestDOY))
      meanSoilTempGas   <- mean(SoilTemperature(seqDescretizationBiogas,localDOYTime,annAvgSurfTemperature,annAvgSurfTempAmplitude,(2*pi/(365)),soilDampingDepthYr,coldestDOY))
      meanSoilTempFloor <- SoilTemperature(digesterHeight,localDOYTime,annAvgSurfTemperature,annAvgSurfTempAmplitude,(2*pi/(365)),soilDampingDepthYr,coldestDOY)
      
      # Calculate the individual heat transfer rates
      qIRRa4 <- viewFactor.cov_sub*QtIntegrated(currentDOY,timeZone,longitude,hr,latitude,P0,masl,aCoef,Sp0,(timeStep/3600))*areaCover*abCover
      qRADa4 <- RadiativeHeatTransfer(sigma,effSkyTemperature,tempSubstrate0,areaCover,surfaceAreaSubstrate,viewFactor.cov_sub,emisCover,emisSubstrate)
      qCON42gr <- surfaceAreaWallBiogas*uOverallBiogas*(meanSoilTempGas - tempSubstrate0)
      qCON4gr <- surfaceAreaWallSubstrate*uOverallSubstrate*(meanSoilTempSub - tempSubstrate0)
      qCONfloorgr <- surfaceAreaWallFloor*uOverallFloor*(meanSoilTempFloor - tempSubstrate0)
      qCONa4 <- areaCover*uOverallCover*(tempAir - tempSubstrate0);
      qADVload4 <- substrateLoadingRate*cpSubstrate*(tempLoading - tempSubstrate0)
      qHeating <- heatingRate
      

      # Compute new substrate
      qTotalSubstrate <- (qADVload4/volumeSubstrate) + (qCONa4/volumeSubstrate) + (qRADa4/volumeSubstrate) + (qIRRa4/volumeSubstrate) + (qCON4gr/volumeSubstrate) + (qCON42gr/volumeSubstrate) + (qCONfloorgr/volumeSubstrate) + (qHeating/volumeSubstrate)
      tempSubstrate <- tempSubstrate0 + ((qTotalSubstrate*timeStep)/(rhoSubstrate*cpSubstrate))
      
      # Insert results into results matrix
      #message("Iterator: ",i)
      mat.res[i,1] <- dos                                             # Insert the Day of Simulation
      mat.res[i,2] <- hr                                              # Insert the hour of the day
      mat.res[i,3] <- localDOYTime                                    # Insert Local DOY Time
      mat.res[i,4] <- tempAir                                         # Insert Temperature of the air                 (K)
      mat.res[i,5] <- tempSubstrate                                   # Insert temperature of the substrate           (K)  
      
      mat.res[i,6] <- qIRRa4
      mat.res[i,7] <- qRADa4
      mat.res[i,8] <- qCON42gr
      mat.res[i,9] <- qCON4gr
      mat.res[i,10] <- qCONfloorgr
      mat.res[i,11] <- qCONa4
      mat.res[i,12] <- qADVload4
      
      mat.res[i,13] <- meanSoilTempSub
      mat.res[i,14] <- meanSoilTempGas
      mat.res[i,15] <- meanSoilTempFloor

      
      # Update to new timestep, temperature and loop contro
      i <- i + 1
      tempSubstrate0 <- tempSubstrate
      hr <- hr + (timeStep/3600)
      
    }
      
  }

  return(mat.res)

}
