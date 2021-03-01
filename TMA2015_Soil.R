# Thermal Model for Anaerobic Digesters (Soil)
# R_TMAD v.0.1
# Written in R version 3.4.0
# Last update: 15/05/2017
# Simon V. Pedersen (May 2017)

####################################################################################################################
#                                                       Functions
####################################################################################################################

##################################################################
# Function: DampingDepth
# Purpose: Calculates the soil-specific damping depth 
# Prog. Id.: SVIP
# Date: 15/05/2017
# Returns: Soil damping depth   (m)
##################################################################

DampingDepth <- function(
  angularFrequency,                        # Angular frequency (i.e. 2 pi/(365*24))         (s^-1)
  soilThermalDiffusivity                   # Thermal Diffusivity of the soil                (m^2 s^-1)
){
  
  # Campell (1977). An introduction to Environmental Biophysics
  dampingDepth <- sqrt((2*soilThermalDiffusivity)/(angularFrequency))
  
  return(dampingDepth)
  
}

##################################################################
# Function: SoilTemperature
# Purpose: Calculates the soil temperature at depth z and time t.
# Prog. Id.: SVIP
# Date: 16/05/2017
# Returns: Soil temperature at depth z and timt t       (k)
##################################################################

SoilTemperature <- function(
  zDepth,                                   # The depth at which soil temperature is calculated                          (m)
  tTime,                                    # The time at which soil temperature is calculated, total.                   (Decimal DOY)
  annAvgSurfTemperature,                    # Annual average soil surface temperature                                    (k)
  annAvgSurfTempAmplitude,                  # Annual average surface temperature fluctuation amplitude                   (k)
  angularFrequency,                         # Angular frequency of the sinusoidal function                               (d^-1)
  dampingDepthYearly,                       # Damping depth                                                              (m)
  coldestDOY                                # The coldest day of year                                                    (Integer)
){
  
  # Campell (1998). An introduction to Environmental Biophysics
  soilTemperature <- annAvgSurfTemperature + ((annAvgSurfTempAmplitude*(sin((angularFrequency*(tTime - coldestDOY))-(pi/2)-(zDepth/dampingDepthYearly))))/(exp(zDepth/dampingDepthYearly)))
  
  return(soilTemperature)
  
}