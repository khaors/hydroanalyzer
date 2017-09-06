#' @title
#' water_budget_direct
#' @description
#' Function to calculate the water budget using the direct method
#' @param Prec Numeric vector with the precipitation values
#' @param EVT Numeric vector with the evapotranspiration values
#' @param Rmax Maximum soil retention (in mm)
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item R: soil storage
#' \item RV: Soil storage variation
#' \item ETR: Real evapotranspiration
#' \item WD: Water deficit
#' \item WE: Water excess
#' \item RN: Runoff
#' \item PC: Deep percolation
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
water_budget_direct <- function(Prec, EVT, Rmax){
  R0 <- Rmax
  ndat <- length(Prec)
  R <- vector('numeric', length = ndat)
  VR <- vector('numeric', length = ndat)
  ETR <- vector('numeric', length = ndat)
  WD <- vector('numeric', length = ndat)
  WE <- vector('numeric', length = ndat)
  RN <- vector('numeric', length = ndat)
  PC <- vector('numeric', length = ndat)
  Rtest <- (Prec[1]-EVT[1])
  for(i in 1:ndat){
    if(i > 1){
      Rtest <- R[i-1] + (Prec[i]-EVT[i])
    }
    else{
      Rtest <- (Prec[i]-EVT[i])
    }
    #print(i)
    if(Rtest > 0 & Rtest <= Rmax){
      R[i] <- Rtest
    }
    else if(Rtest > Rmax){
      R[i] <- Rmax
    }
    else if(Rtest < 0){
      R[i] <- 0
    }
    #print(c(Prec[i]-EVT[i],Rtest,R[i]))
    # Variation in soil storage
    if(i == 1){
      VR[i] <- R[i]
    }
    else {
      VR[i] <- R[i]-R[i-1]
    }
    # Real evapotranspiration
    water_season <- Prec[i]-EVT[i]
    if(water_season > 0){
      ETR[i] <- EVT[i]
    }
    else if(water_season <= 0){
      ETR[i]<- Prec[i]+abs(VR[i])
    }
    # Water deficit
    WD[i] <- EVT[i]-ETR[i]
    # Water Excess
    if(water_season >0){
      WE[i] <- water_season - VR[i]
    }
    else{
      WE[i] <- 0
    }
    # Runoff and Deep Percolation
    if(i == 1){
      RN[i] <- 0.5*WE[i]
      PC[i] <- 0.5*WE[i]
    }
    else {
      RN[i]<-0.5*(RN[i-1]+WE[i])
      PC[i] <- WE[i]-RN[i]
    }
  }
  results <- list(R = R, VR = VR, WD = WD, WE = WE, RN = RN, PC = PC, ETR = ETR)
}
#' @title
#' evt_thornthwaite
#' @description
#' Function to calculate the potential evapotranspiration using the Thornthwaite's approach
#' @param Temp A numeric vector with the mean monthly temperature
#' @return
#' A numeric vector with the calculate EVT
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
evt_thornthwaite <- function(Temp){
  ndat <- length(Temp)
  if(ndat != 12){
    stop('ERROR: The number of montly averages is not equal to 12')
  }
  I_month <- vector('numeric', length = ndat)
  I_month <- (Temp/5)**1.514
  I_year <- sum(I_month)
  alpha <-  6.75e-7*I_year**3-7.71e-5*I_year**2+0.0179*I_year+0.49
  evt <- 1.6*(10*Temp/I_year)**alpha
  return(evt)
}
