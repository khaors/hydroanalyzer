#' @title
#' abcd.year.model
#' @description
#' Function to calculate the components of the hydrological cycle using the abcd model on a yearly basis
#' @param par.model A numeric vector of length 5 with the values of the abcd parameters and the initial groundwater storage
#' @param Prec A numeric vector with the values of the Precipitation (on a yearly basis)
#' @param PEV A numeric vector with the values of potential evapotranspiration
#' @param ... Additional parameters (used for compatibility reasons)
#' @return
#' A list with the following elements:
#' \itemize{
#' \item Qt: Numeric vector with the yearly discharges
#' \item SRt: Numeric vector with the surface runnoff
#' \item ETt: Numeric vector with the actual evapotranspiration
#' \item Dt: Numeric vector with the deep percolation
#' \item BFt: Numeric vector with the base flow
#' \item GFt: Numeric vector with the groundwater flow in the aquifer
#' \item GSt: Numeric vector with the groundwater storage
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family abcd_model functions
#' @export
abcd.year.model <- function(par.model = c(0,0,0,0,0), Prec, PEV, ...){
  #
  a <- par.model[1]
  b <- par.model[2]
  c <- par.model[3]
  d <- par.model[4]
  # Initial GW storage
  G_1 <- par.model[5]
  # Initialize output vectors
  number_years <- length(Prec)
  GS <- vector('numeric', length = number_years)
  BF <- vector('numeric', length = number_years)
  Qs <- vector('numeric', length = number_years)
  GF <- vector('numeric', length = number_years)
  #
  #Surface runnof
  SR <- a*Prec
  In <- (1.0-a)*Prec
  #Evapotranspiration
  E <- b*In
  # Deep percolation
  DP <- (1.0-b)*In
  #Base flow
  BF[1] <- c*G_1
  # GW Flow
  GF[1] <- d*G_1
  for(it in 1:length(Prec)){
    if(it == 1){
      BF[it] <- c*G_1
      GF[it] <- d*G_1
      GS[it] <- (G_1-BF[it]-GF[it])+DP[it]
    }
    else {
      BF[it] <- c*GS[it-1]
      GF[it] <- d*GS[it-1]
      GS[it] <- (GS[it-1]-BF[it]-GF[it])+DP[it]
    }
    Qs[it] <- SR[it] + BF[it]
  }
  #res <- Qs
  res <- list(Qt = Qs, SRt = SR, ETt = E, Dt = DP, BFt = BF, GSt = GS)
  return(res)
}
#' @title
#' abcd.month.model
#' @description
#' Function to calculate the components of the hydrological cycle using the abcd model on a monthly basis
#' @param par.model A numeric vector of length 5 with the values of the abcd parameters and the initial groundwater storage
#' @param Prec A numeric vector with the values of the Precipitation (on a yearly basis)
#' @param PEV A numeric vector with the values of potential evapotranspiration
#' @param ... Additional parameters (used for compatibility reasons)
#' @return
#' A list with the following elements:
#' \itemize{
#' \item Qt: Numeric vector with the yearly discharges
#' \item SRt: Numeric vector with the surface runnoff
#' \item Wt: Numeric vector witht the available soil water
#' \item SMt: Numeric vector with the soil moisture
#' \item ETt: Numeric vector with the actual evapotranspiration
#' \item Dt: Numeric vector with the deep percolation
#' \item BFt: Numeric vector with the base flow
#' \item GFt: Numeric vector with the groundwater flow in the aquifer
#' \item GSt: Numeric vector with the groundwater storage
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family abcd_model functions
#' @export
abcd.month.model <- function(par.model = c(0,0,0,0,0,0,0), Prec, PEV, ...){
  a <- par.model[1]
  b <- par.model[2]
  c <- par.model[3]
  d <- par.model[4]
  #
  G_1 <- par.model[5]
  S_1 <- par.model[6]
  #FC <- par.model[7]
  #
  Ep <- PEV
  number_months <- length(Prec)
  Qt <- vector('numeric', length = number_months)
  Wt <- vector('numeric', length = number_months)
  St <- vector('numeric', length = number_months)
  Rt <- vector('numeric', length = number_months)
  Gt <- vector('numeric', length = number_months)
  Qbt <- vector('numeric', length = number_months)
  #
  Q <- 0.0
  for(irow in 1:number_months){
    #Available soil water
    W <- Prec[irow]+S_1
    Wt[irow] <- W
    #Evapotranspiration potential
    yy <- (W+b)/(2*a)-(((W+b)/(2*a))^2-((W*b)/a))^0.5
    #Potential evapotranspiration
    E <- yy*(1-exp(-Ep[irow]/b))
    #Soil moisture
    S <- yy- E
    St[irow] <- S
    #Runoff
    Qd <- (1-c)*(W-yy)
    #GW Recharge
    R <- c*(W-yy)
    Rt[irow] <- R
    # Groundwater storage
    G <- (1/(1+d))*(R+G_1)
    Gt[irow] <- G
    # Base flow
    Qb <- d*G
    Qbt[irow] <- Qb
    # Total discharge
    Q <- Qb+Qd
    G_1 <- G
    S_1 <- S
    Qt[irow] <- Q
  }
  res <- list(Qt = Qt, Wt = Wt, SRt = Qd, SMt = St, Rt = Rt, Gt = Gt,
              Qbt = Qbt, ETt = E)
  #res <- Qt
  return(res)
}
#' @title
#' abcd.mod.year.model
#' @description
#' Function to calculate the components of the hydrological cycle using the modified abcd model
#' on a yearly basis. The modification of the original abcd model takes into account the pumping
#' rates and the depth to the aquifer used to calculate the variations in the groundwater levels
#' @param par.model A numeric vector of length 5 with the values of the abcd parameters and the initial groundwater storage
#' @param Prec A numeric vector with the values of the Precipitation (on a yearly basis)
#' @param PEV A numeric vector with the values of potential evapotranspiration
#' @param Vp A numeric vector with the values of the volumes of water extracted by pumping
#' @param ... Additional parameters (used for compatibility reasons)
#' @return
#' A list with the following elements:
#' \itemize{
#' \item Qt: Numeric vector with the yearly discharges
#' \item SRt: Numeric vector with the surface runnoff
#' \item ETt: Numeric vector with the actual evapotranspiration
#' \item Dt: Numeric vector with the deep percolation
#' \item BFt: Numeric vector with the base flow
#' \item GFt: Numeric vector with the groundwater flow in the aquifer
#' \item GSt: Numeric vector with the groundwater storage
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family abcd_model functions
#' @export
abcd.mod.year.model <- function(par.model = c(0,0,0,0,0), Prec, PEV, Vp, ...){
  #
  a <- par.model[1]
  b <- par.model[2]
  c <- par.model[3]
  d <- par.model[4]
  # Initial GW storage
  G_1 <- par.model[5]
  # Initialize output vectors
  number_years <- length(Prec)
  GS <- vector('numeric', length = number_years)
  BF <- vector('numeric', length = number_years)
  Qs <- vector('numeric', length = number_years)
  GF <- vector('numeric', length = number_years)
  Vp1 <- NULL
  if(length(Vp) == 1){
    Vp1 <- vector('numeric', length = number_years)
    Vp1[1:number_years] <- Vp
  }
  else {
    Vp1 <- Vp
  }
  #
  #Surface runnof
  SR <- a*Prec
  In <- (1.0-a)*Prec
  #Evapotranspiration
  E <- b*In
  # Deep percolation
  DP <- (1.0-b)*In
  #Base flow
  BF[1] <- c*G_1
  # GW Flow
  GF[1] <- d*G_1
  for(it in 1:length(Prec)){
    if(it == 1){
      BF[it] <- c*G_1
      GF[it] <- d*G_1
      GS[it] <- (G_1-BF[it]-GF[it]-Vp1[it])+DP[it]
    }
    else {
      BF[it] <- c*GS[it-1]
      GF[it] <- d*GS[it-1]
      GS[it] <- (GS[it-1]-BF[it]-GF[it]-Vp1[it])+DP[it]
    }
    Qs[it] <- SR[it] + BF[it]
  }
  #res <- Qs
  res <- list(Qt = Qs, SRt = SR, ETt = E, Dt = DP, BFt = BF, GSt = GS)
  return(res)
}
#' @title
#' abcd.mod.month.model
#' @description
#' Function to calculate the components of the hydrological cycle using the modified abcd model
#' on a monthly basis. The modification of the original abcd model takes into account the pumping
#' rates and the depth to the aquifer used to calculate the variations in the groundwater levels
#' @param par.model A numeric vector of length 5 with the values of the abcd parameters and the initial groundwater storage
#' @param Prec A numeric vector with the values of the Precipitation (on a yearly basis)
#' @param PEV A numeric vector with the values of potential evapotranspiration
#' @param Vp A numeric vector with the values of the volumes of water extracted by pumping
#' @param ... Additional parameters (used for compatibility reasons)
#' @return
#' A list with the following elements:
#' \itemize{
#' \item Qt: Numeric vector with the yearly discharges
#' \item SRt: Numeric vector with the surface runnoff
#' \item Wt: Numeric vector witht the available soil water
#' \item SMt: Numeric vector with the soil moisture
#' \item ETt: Numeric vector with the actual evapotranspiration
#' \item Dt: Numeric vector with the deep percolation
#' \item BFt: Numeric vector with the base flow
#' \item GFt: Numeric vector with the groundwater flow in the aquifer
#' \item GSt: Numeric vector with the groundwater storage
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family abcd_model functions
#' @export
abcd.mod.month.model <- function(par.model = c(0,0,0,0,0,0,0), Prec, PEV, Vp, ...){
  a <- par.model[1]
  b <- par.model[2]
  c <- par.model[3]
  d <- par.model[4]
  #
  G_1 <- par.model[5]
  S_1 <- par.model[6]
  #FC <- par.model[7]
  #
  Ep <- PEV
  number_months <- length(Prec)
  Qt <- vector('numeric', length = number_months)
  Wt <- vector('numeric', length = number_months)
  St <- vector('numeric', length = number_months)
  Rt <- vector('numeric', length = number_months)
  Gt <- vector('numeric', length = number_months)
  Qbt <- vector('numeric', length = number_months)
  Vp1 <- NULL
  if(length(Vp) == 1){
    Vp1 <- vector('numeric', length = number_months)
    Vp1[1:number_months] <- Vp
  }
  else {
    Vp1 <- Vp
  }
  #
  Q <- 0.0
  for(irow in 1:number_months){
    #Available soil water
    W <- Prec[irow]+S_1
    Wt[irow] <- W
    #Evapotranspiration potential
    yy <- (W+b)/(2*a)-(((W+b)/(2*a))^2-((W*b)/a))^0.5
    #Potential evapotranspiration
    E <- yy*(1-exp(-Ep[irow]/b))
    #Soil moisture
    S <- yy- E
    St[irow] <- S
    #Runoff
    Qd <- (1-c)*(W-yy)
    #GW Recharge
    R <- c*(W-yy)
    Rt[irow] <- R
    # Groundwater storage
    G <- (1/(1+d))*(R+G_1)
    Gt[irow] <- G-Vp1[irow]
    # Base flow
    Qb <- d*G
    Qbt[irow] <- Qb
    # Total discharge
    Q <- Qb+Qd
    G_1 <- G
    S_1 <- S
    Qt[irow] <- Q
  }
  res <- list(Qt = Qt, Wt = Wt, SRt = Qd, SMt = St, Rt = Rt, Gt = Gt,
              Qbt = Qbt, ETt = E)
  #res <- Qt
  return(res)
}
