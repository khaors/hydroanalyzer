#' @title
#' nathan_filter
#' @description
#' Function to calculate the baseflow from a discharge time series using the digital filter
#' proposed by Nathan1991.
#' @param Q A numeric vector with the discharge time series
#' @param a A numeric value specifying the filter parameter
#' @return
#' A numeric vector with the calculated baseflow
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family digital_filter functions
#' @export
#' @source Nathan, R. and McMahon, T. 1991. Evaluation of automated techniques for base flow
#' and recession analyses, Water Resources Research, 27, 1783-1784
nathan_filter <- function(Q, a){
  nq <- length(Q)
  Qbf <- vector('numeric', length = nq)
  Qbf[1] <- 0.5*Q[1]
  for(it in 2:nq){
    Qbf[it] <- a*Q[it-1]+0.5*(1-a)*(Q[it]+Q[it-1])
    if(Qbf[it] > Q[it]){
      Qbf[it] <- Q[it]
    }
  }
  return(Qbf)
}
#' @title
#' chapman_filter
#' @description
#' Function to calculate the baseflow from a discharge time series using the digital filter
#' proposed by Chapman1990.
#' @param Q A numeric vector with the discharge time series
#' @param a A numeric value specifying the filter parameter
#' @return
#' A numeric vector with the calculated baseflow
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family digital_filter functions
#' @export
#' @source Chapman, T. G. 1991. Comment on "Evaluation of automated techniques for base flow
#' and recession analyses" by R. J. Nathan and T. A. McMahon, Water Resources Research, 27,
#' 1783-1784
chapman_filter <- function(Q, a){
  nq <- length(Q)
  Qbf <- vector('numeric', length = nq)
  Qbf[1] <- 0.5*Q[1]
  coef1 <- (3*a-1)/(3-a)
  coef2 <- (1-a)/(3-a)
  for(it in 2:nq){
    Qbf[it] <- coef1*Qbf[it-1] + coef2*(Q[it]+Q[it-1])
    #Qbf[it] <- (a/(2-a))*Qbf[it-1]+((1-a)/(2-a))*Q[it]
    if(Qbf[it] > Q[it]){
      Qbf[it] <- Q[it]
    }
  }
  return(Qbf)
}
#' @title
#' eckhardt_filter
#' @description
#' Function to calculate the baseflow from a discharge time series using the digital filter
#' proposed by Eckhardt2005.
#' @param Q A numeric vector with the discharge time series
#' @param a A numeric value specifying the filter parameter
#' @param BFImax A numeric value specifying the maximum Base Flow Index used in the separation.
#' @return
#' A numeric vector with the calculated baseflow
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family digital_filter functions
#' @export
#' @source Eckhardt, K.,  2005. How to construct recursive digital filters for baseflow
#' separation, Hydrological Processes, 19, 507-515
eckhardt_filter <- function(Q, a, BFImax){
  nq <- length(Q)
  Qbf <- vector('numeric', length = nq)
  Qbf[1] <- 0.5*Q[1]
  coef1 <- (1-BFImax)*a
  coef2 <- (1-a)*BFImax
  coef3 <- (1-a*BFImax)
  for(it in 2:nq){
    Qbf[it] <- (coef1*Qbf[it-1]+coef2*Q[it])/coef3
    if(Qbf[it] > Q[it]){
      Qbf[it] <- Q[it]
    }
  }
  return(Qbf)
}
