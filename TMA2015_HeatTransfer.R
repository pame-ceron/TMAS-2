# Thermal Model for Anaerobic Digesters (HeatTransfer)
# R_TMAD v.0.1
# Written in R version 3.4.0
# Last update: 15/05/2017
# Simon V. Pedersen (May 2017)

####################################################################################################################
#                                                       Functions
####################################################################################################################

##################################################################
# Function: OverallHeatTransferCoef
# Purpose: Calculates the overall heat transfer coefficient 
# Prog. Id.: SVIP
# Date: 15/05/2017
# Returns: The overll heat transfer coefficient (W m^-2 K^-1)
##################################################################

OverallHeatTransferCoef <- function(
  frameResistance                       # The data.frame containing the resistances as (deltaL,k) or (1,h)
){
  
  # Calculate the value of the resistance, as thickness/(h or k)
  list.Resistances <- mapply(function(thickness,resist) {thickness/resist},frameResistance[[1]],frameResistance[[2]])
  
  # Calculate the overall heat transfer coefficient U, as the reciprocal of the sum of list.resistances.
  overallU <- 1/(sum(list.Resistances))
  
  return(overallU)
  
}

##################################################################
# Function: RadiativeHeatTransfer
# Purpose: Calculates the overall heat transfer coefficient 
# Prog. Id.: SVIP
# Date: 15/05/2017
# Returns: The overll heat transfer coefficient (W m^-2 K^-1)
##################################################################

RadiativeHeatTransfer <- function(
  stefanBoltzmannConstant,          # The Stefan-Boltzmann constant, sigma                  (W m^-2 K^-4)
  skyTemperature,                   # The effective sky temperature                         (K)
  slurryTemperature,                # The slurry temperature                                (K)
  areaCover,                        # Surface-area of the cover exposed to the atmosphere   (m^2)
  surfaceAreaSlurry,                # Surface area of the gas-slurry interface              (m^2)
  viewFactor.cov.sub,               # View factor of from cover to substrate (slurry)       (DIMENSIONLESS)
  emisCover,                        # Emissivity of the cover                               (DIMENSIONLESS)
  emisSlurry                        # EMissivity of the substrate (slurry)                  (DIMENSIONLESS)
){
  
  # Equation simplified from the one given by Cengel (2007). Heat and Mass Transfer: A practical approach. 
  # Assumptions:
  #   + The sky is a perfect blackbody, and its emissivity = 1.
  #   + The over- and underside of the cover is of the same radiative characteristics. I.e. same emissivity.
  #   + Applicability of the reciprocity relation, to avoid A_sky.
  
  numerator <- stefanBoltzmannConstant*((skyTemperature^(4))-(slurryTemperature^(4)))
  denominator <- (1/areaCover) + (2*((1-emisCover)/(areaCover*emisCover))) + (1/(areaCover*viewFactor.cov.sub)) + ((1-emisSlurry)/(surfaceAreaSlurry*emisSlurry))
  
  qRadiative.sky.sub <- (numerator)/(denominator)
  
  return(qRadiative.sky.sub)

}