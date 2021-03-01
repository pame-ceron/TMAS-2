# Thermal Model for Anaerobic Digesters (Geometry)
# R_TMAD v.0.1
# Written in R version 3.4.0
# Last update: 17/05/2017
# Simon V. Pedersen (May 2017)

####################################################################################################################
#                                                       Functions
####################################################################################################################

##################################################################
# Function: CylindricalViewFactor
# Purpose: Calculates the view factor for two coaxial parallel 
#          disks.
# Prog. Id.: SVIP
# Date: 17/05/2017
# Returns: The view factor Fij    (DIMENSIONLESS)
##################################################################

CylendricalViewFactor <- function(
  radiusCover,                            # The radius of the cover                            (m)
  radiusSubstrateInterface,               # The radius of the substrate-biogas interface       (m)
  distanceSeparation                      # Distance of separation between the two disks       (m) 
){
  
  # The following calculations are obtained from: Cengel (2007). Heat and Mass Transfer: A practical approach.
  ri <- (radiusCover/distanceSeparation)
  rj <- (radiusSubstrateInterface/distanceSeparation)
  sVal <- 1 + ((1+(rj^2))/(ri^2))
  fij <- (1/2)*(sVal-(((sVal^2)-(4*((rj/ri)^2)))^(1/2)))
  
  return(fij)
  
}















