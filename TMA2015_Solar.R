# Thermal Model for Anaerobic Digesters (Solar)
# R_TMAD v.0.1
# Written in R version 3.4.0
# Last update: 15/05/2017
# Simon V. Pedersen (May 2017)

####################################################################################################################
#                                                       Functions
####################################################################################################################

##################################################################
# Function: Theta
# Purpose: Calculating the suns elevation angle
# Prog. Id.: SVIP
# Date: 10/05/2017
# Returns: The Sun's elevation angle (DEGREE)
##################################################################

Theta <- function(
  latitude,                       # Location Latitude      (Decimal Degrees)
  declinationAngle,               # Declination Angle      (Decimal Degrees)
  hourAngle                       # Hour Angle             (Decimal Degrees)
){
  
  # Degrees to Radian conversion
  deg.rad <- (pi/180)
  # Radian to Degree conversion
  rad.deg <- (180/pi)

  # Campbell (1998). An Introduction to Envuronmental Biophysics  
  elevAngle <- asin((sin(latitude*deg.rad)*sin(declinationAngle*deg.rad)) + (cos(latitude*deg.rad)*cos(declinationAngle*deg.rad)*cos(hourAngle*deg.rad)))*rad.deg

  return(elevAngle)
  
}

##################################################################
# Function: DeclinationAngle
# Purpose: Calculating the suns' declination angle
# Prog. Id.: SVIP
# Date: 10/05/2017
# Returns: The sun's declination angle (DEGREE)
##################################################################

DeclinationAngle <- function(   
  doy                             # The Day of Year (Integer)
){
  
  # Cooper (1969). The Absorption of Radiation in Solar Stills.
  decAngle <- 23.45*sin(((2*pi*(284+doy))/(365)))
  
  return(decAngle)
  
}

##################################################################
# Function: AirMass
# Purpose: Calculates the air-mass number. A relative measure of
#          distance traveled by the incident beam through the
#          atmosphere.
# Prog. Id.: SVIP
# Date: 10/05/2017
# Returns: The Optical Air Mass Number (Dimensionless)
##################################################################

AirMass <- function(   
  pa,                             # Pressure at the digester location                 (Pa)
  p0,                             # Standard atmosphere. Pressure at sea-level        (Pa)
  elevationAngle                  # The sun's elevation angle                         (RADIANS)
){
  
  # Degrees to Radian conversion
  deg.rad <- (pi/180)
  
  # Gebremedhin et al. (2005). Heat Transfer Model for Plug Flow Anaerobic Digesters.
  optAirMass <- ((pa/p0)/(sin(elevationAngle*deg.rad)))
  
  return(optAirMass)
  
}
  
##################################################################
# Function: Sb
# Purpose: Calculates direct irradiance on a horizontal surface
# Prog. Id.: SVIP
# Date: 10/05/2017
# Returns: Direct irradiance on horizontal surface (W m^-2)
##################################################################

Sb <- function(   
  elevationAngle,                 # The sun's elevation angle                                                   (DEGREE)
  sp                              # Direct solar irradiance on a surface perpendicular to the solar beam        (W m^-2)
){
  
  # Degrees to Radian conversion
  deg.rad <- (pi/180)
  
  # Campbell (1977). An Introduction to Environmental Biophysics
  sbout <- sp*sin(elevationAngle*deg.rad)
  
  return(sbout)
}

##################################################################
# Function: Sp
# Purpose: Calculates the direct solar irradiance on a surface
#          perendicular to the solar beam.
# Prog. Id.: SVIP
# Date: 10/05/2017
# Returns: Direct irradiance on perpendicular surface (W m^-2)
##################################################################

Sp <- function(   
  aCoef,                          # Coefficient of transmissivity. Constant                         (Dimensionless)
  sp0,                            # Extraterrestrial flux density normal to the solar beam          (W m^-2)
  airMassNo                       # Optical Air-mass number                                         (Dimensionless)
){
  
  # Campbell (1977). An Introduction to Environmental Biophysics
  spout <- (aCoef^(airMassNo))*sp0
  
  return(spout)
  
}  
  
##################################################################
# Function: Sd
# Purpose: Calculates diffuse sky irradiance
# Prog. Id.: SVIP
# Date: 10/05/2017
# Returns: Diffuse sky irradiance   (W m^-2) 
##################################################################

Sd <- function(   
  aCoef,                          # Coefficient of transmissivity. Constant                         (Dimensionless)
  sp0,                            # Extraterrestrial flux density normal to the solar beam          (W m^-2)
  airMassNo,                      # Optical Air-mass number                                         (Dimensionless)
  elevationAngle                  # The suns' elevation angle                                       (DEGREES)
){
  
  # Degrees to Radian conversion
  deg.rad <- (pi/180)
  
  # Campbell and Norman (1998). An Introduction to Environmental Biophysics
  sdout <- 0.3*sp0*(1-(aCoef^(airMassNo)))*sin(elevationAngle*deg.rad)
  
  return(sdout)
  
}

##################################################################
# Function: TimeCorrectionFactor
# Purpose: Calculates the time correction factor
# Prog. Id.: SVIP
# Date: 15/05/2017
# Returns: The time correction factor.   (min)
# NOTE: The following calculation of the time correction factor is 
#       based on the method outlined by the photovoltaic education
#       education network, at www.pveducation.org/pvcdrom. The
#       following function has been re-verified with PVEducation
#       on the 10th May 2017 by SVIP.
##################################################################

TimeCorrectionFactor <- function(   
  timeZone,                             # Timezone fromt GMT               (hr)
  doy,                                  # Day of year                      (Integer)
  longitude                            # Location longitude               (DEGREES)
){
  
  # Degrees to Radian conversion
  deg.rad <- (pi/180)
  # Radian to Degree conversion
  rad.deg <- (180/pi)
  
  lstm <- 15*(timeZone)                                             # Local Standard Time Meridian                   (DEGREES)
  bVal <- ((2*pi)/(365))*(doy-81)                                   # B-value, used for the Equation of Time         (RADIAN)
  eot <- (9.87*sin(2*bVal))-(7.53*cos(bVal))-(1.5*sin(bVal))        # The Equation of Time                           (min)
  
  # The time Correction Factor considers the difference in Local Solar Time, due to longitudinal variation within the time-zone
  tcFactor <- (4*(longitude-lstm))+eot                              # The Time Correction Factor                     (min)
  
  return(tcFactor)
  
}

##################################################################
# Function: HourAngle
# Purpose: Calculates the Hour Angle of the sun
# Prog. Id.: SVIP
# Date: 10/05/2017
# Returns: The Hour Angle of the sun.   (DEGREES)
# NOTE: The following calculation of the sun's hour angle, is 
#       based on the method outlined by the photovoltaic education
#       education network, at www.pveducation.org/pvcdrom. The
#       following function has been re-verified with PVEducation
#       on the 10th May 2017 by SVIP.
##################################################################

HourAngle <- function(   
  timeCorrection,                       # The Time Correction Factor       (min)
  localTime                             # Local time at the location       (hr)
){
  
  # Degrees to Radian conversion
  deg.rad <- (pi/180)
  # Radian to Degree conversion
  rad.deg <- (180/pi)
 
  lst <- localTime + (timeCorrection/60)                            # The Local Solar Time                           (hr)
  hourAngleOut <- 15*(lst-12)                                       # The Hour Angle                                 (DEGREE)
  
  return(hourAngleOut)
  
}

##################################################################
# Function: TimeSunrise
# Purpose: Calculates the local time of sunrise
# Prog. Id.: SVIP
# Date: 11/05/2017
# Returns: The local time of sunrise   (hr)
##################################################################

TimeSunrise <- function(   
  latitude,                             # Location latitude              (DEGREES)
  declinationAngle,                     # The suns declination angle     (DEGREES)
  tcFactor                              # The Time Correction Factor     (min)
){
  
  # Degrees to Radian conversion
  deg.rad <- (pi/180)
  # Radian to Degree conversion
  rad.deg <- (180/pi)
  
  # Based on information from The Photovoltaic Education Network, PVEducation.org. Verified 11th May 2017. 
  tSunrise <- 12 - (((1)/(15*deg.rad))*acos(((-sin(latitude*deg.rad)*sin(declinationAngle*deg.rad))/(cos(latitude*deg.rad)*cos(declinationAngle*deg.rad))))) - (tcFactor/60)
  
  return(tSunrise)
  
}

##################################################################
# Function: TimeSunset
# Purpose: Calculates the local time of sunset
# Prog. Id.: SVIP
# Date: 11/05/2017
# Returns: The local time of sunset   (hr)
##################################################################

TimeSunset <- function(   
  latitude,                             # Location latitude              (DEGREES)
  declinationAngle,                     # The suns declination angle     (DEGREES)
  tcFactor                              # The Time Correction Factor     (min)
){
  
  # Degrees to Radian conversion
  deg.rad <- (pi/180)
  # Radian to Degree conversion
  rad.deg <- (180/pi)
  
  # Based on information from The Photovoltaic Education Network, PVEducation.org. Verified 11th May 2017. 
  tSunset <- 12 + (((1)/(15*deg.rad))*acos(((-sin(latitude*deg.rad)*sin(declinationAngle*deg.rad))/(cos(latitude*deg.rad)*cos(declinationAngle*deg.rad))))) - (tcFactor/60)
  
  return(tSunset)
  
}

##################################################################
# Function: PressureAtAltitude
# Purpose: Calculates the atmospheric pressure at m.a.s.l.
# Prog. Id.: SVIP
# Date: 11/05/2017
# Returns: Atmospheric pressure at m.a.s.l.   (Pa)
##################################################################

PressureAtAltitude <- function(   
  p0,                                   # Standard atmosphere at 0 m.a.sl.              (Pa)
  masl                                  # Meters above sea-level                        (m)
){
  
  pAltitude <- p0*exp(masl/8000)
  
  return(pAltitude)
  
}

##################################################################
# Function: Qt
# Purpose: Calculates the total irradiance on a horizontal surface
# Prog. Id.: SVIP
# Date: 11/05/2017
# Returns: Total irradiance on horizontal surface   (W m^-2)
##################################################################

Qt <- function(   
  doy,                                  # Day of year                                                                  (Integer)
  timeZone,                             # Location time-zone from GMT                                                  (hr)
  longitude,                            # Location longitude                                                           (DEGREES)
  localTime,                            # Local time at the location                                                   (hr)
  latitude,                             # Location latitude                                                            (DEGREES)
  p0,                                   # Atmospheric pressure at 0 m.a.s.l                                            (Pa)
  masl,                                 # Digester location above sea-level                                            (m)
  aCoef,                                # Coefficient of Tranismissivity. Constant. Average annual condition           (Dimensionless)
  sp0                                   # Extraterrestrial flux density normal to the solar beam                       (W m^-2)
){
  
  # Degrees to Radian conversion
  deg.rad <- (pi/180)
  # Radian to Degree conversion
  rad.deg <- (180/pi)
  
  # Calculate Intermediate variables
  valDeclination <- DeclinationAngle(doy)
  valTimeCorrection <- TimeCorrectionFactor(timeZone,doy,longitude)
  valHourAngle <- HourAngle(valTimeCorrection,localTime)
  valTheta <- Theta(latitude,valDeclination,valHourAngle)
  
  if(valTheta < 0){
    valTheta <- 0
  }
  
  valPa <- PressureAtAltitude(p0,masl)
  valAirMass <- AirMass(valPa,p0,valTheta)
  
  #message("FROM QT - valDeclination:", valDeclination)
  #message("FROM QT - valTimeCorrection:", valTimeCorrection)
  #message("FROM QT - valHourAngle:", valHourAngle)
  #message("FROM QT - valTheta:", valTheta)
  #message("FROM QT - valPa:", valPa)
  #message("FROM QT - valAirMass:", valAirMass)
  
  # Based on information from The Photovoltaic Education Network, PVEducation.org. Verified 11th May 2017. 
  qTotal <- Sb(valTheta,Sp(aCoef,sp0,valAirMass))+Sd(aCoef,sp0,valAirMass,valTheta)
  
  #message("FROM QT - Sp:", Sp(aCoef,sp0,valAirMass))
  #message("FROM QT - Sb:", Sb(valTheta,Sp(aCoef,sp0,valAirMass)))
  #message("FROM QT - Sd:", Sd(aCoef,sp0,valAirMass,valTheta))
  
  return(qTotal)
  
}

##################################################################
# Function: QtIntegrated
# Purpose: Calculates the total irradiance on a horizontal surface
#          Not really integrated, but the flux mid-way of two time
#          points.
# Prog. Id.: SVIP
# Date: 15/05/2017
# Returns: Integrated total irradiation on horizontal plane(W m^-2)
##################################################################

QtIntegrated <- function(   
  doy,                                  # Day of year                                                                  (Integer)
  timeZone,                             # Location time-zone from GMT                                                  (hr)
  longitude,                            # Location longitude                                                           (DEGREES)
  localTime,                            # Local time at the location                                                   (hr,decimal)
  latitude,                             # Location latitude                                                            (DEGREES)
  p0,                                   # Atmospheric pressure at 0 m.a.s.l                                            (Pa)
  masl,                                 # Digester location above sea-level                                            (m)
  aCoef,                                # Coefficient of Tranismissivity. Constant. Average annual condition           (Dimensionless)
  sp0,                                  # Extraterrestrial flux density normal to the solar beam                       (W m^-2)
  deltaTimeT                            # The time-step taken                                                          (hr)            
){
  
  # Degrees to Radian conversion
  deg.rad <- (pi/180)
  # Radian to Degree conversion
  rad.deg <- (180/pi)
  
  qTotalIntegrated <- -999
  
  # Calculate Intermediate variables
  valDeclination <- DeclinationAngle(doy)
  valTimeCorrection <- TimeCorrectionFactor(timeZone,doy,longitude)
  probeTimeT <- (localTime - (deltaTimeT/2))        # Evaluate at the time midway between the previous time-point and the current time-point
  timeIntervalLocalTime <- (localTime - deltaTimeT)
  tSunrise <- TimeSunrise(latitude,valDeclination,valTimeCorrection)
  tSunset <- TimeSunset(latitude,valDeclination,valTimeCorrection)
  
  #message("FROM SOLAR, qTotalIntegrated: ", qTotalIntegrated)
  #message("FROM SOLAR, DOY: ", doy)
  #message("FROM SOLAR, timeZone: ", timeZone)
  #message("FROM SOLAR, longitude: ", longitude)
  #message("FROM SOLAR, localTime: ", localTime)
  #message("FROM SOLAR, latitude: ", latitude)
  #message("FROM SOLAR, p0: ", p0)
  #message("FROM SOLAR, masl: ", masl)
  #message("FROM SOLAR, aCoef: ", aCoef)
  #message("FROM SOLAR, sp0: ", sp0)
  #message("FROM SOLAR, deltaTimeT: ", deltaTimeT)

  #message("QTIntegrated: ", timeIntervalLocalTime)
  #message("QTIntegrated: ", tSunrise)
  #message("QTIntegrated: ", tSunset)
  
  if((timeIntervalLocalTime < tSunrise) || (timeIntervalLocalTime > tSunset)){
    qTotalIntegrated <- 0
  }else{
    qTotalIntegrated <- Qt(doy,timeZone,longitude,probeTimeT,latitude,p0,masl,aCoef,sp0)
  }
  
  #message("FROM SOLAR, qTotalIntegrated: ", qTotalIntegrated)
  
  if(qTotalIntegrated < 0){
    qTotalIntegrated <- 0
  }
  
  return(qTotalIntegrated)
  
}



  
  
  
  