# Thermal Model for Anaerobic Digesters (AmbientAtmosphere)
# R_TMAD v.0.1
# Written in R version 3.4.0
# Last update: 15/05/2017
# Simon V. Pedersen (May 2017)

####################################################################################################################
#                                                       Functions
####################################################################################################################

##################################################################
# Function: TSky
# Purpose: Calculates the effective sky temperature 
# Prog. Id.: SVIP
# Date: 15/05/2017
# Returns: The effective sky temperature (Kelvin)
##################################################################

SkyTemperature <- function(
  airTemperature                       # Ambient air temperature (Kelvin)
){
  
  # Derived from Swinbank (1963)
  tSky <- 0.0552*(airTemperature^(3/2))
  
  return(tSky)
  
}

##################################################################
# Function: AmbientAirTemperature
# Purpose: Calculate the ambient air temperature, based on a sine-
#          approximation of the annual temperature fluctiations.
# Prog. Id.: SVIP
# Date: 15/05/2017
# Returns: The ambient air temperature  (K)
##################################################################

AmbientAirTemperature <- function(
  annualAverageTemperature,                      # Annual average air temperature                     (K)
  annualTemperatureAmplitude,                    # Annual amplitude of temperature fluctuations       (K)
  coldestDOY,                                    # Coldest Day of Year                                (Integer)
  localDOYTime                                   # The current DOY and HR as a decimal                (Decimal)
){

  airTemperature <- annualAverageTemperature + (annualTemperatureAmplitude*sin((((2*pi)/(365))*(localDOYTime-coldestDOY))-(pi/2)))
  
  return(airTemperature)
  
}