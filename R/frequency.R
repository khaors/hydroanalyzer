#' @title
#' empirical_frequency
#' @description
#' Function to calculate the empirical frequency associated with a given dataset using a specific formula.
#' @param Var1 A numeric vector with the values of the hydrologic variable to be analized
#' @param model A character string specifying the model used to calculate the empirical frequency
#' @return
#' A list with the following entries
#' \itemize{
#' \item Var: A numeric vector with the sorted values of the hydrologic variable
#' \item Freq: A numeric vector with the cumulative frequency associated with each value of the
#' hydrologic variable
#' \item Treturn: Return period associated with each value of the hydrological variable calculated
#' using the specified model
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
empirical_frequency <- function(Var1, model){
  model.types <- c("None", "weibull", "median", "hosking",
                   "blom", "cunnane",
                   "gringorten", "hazen")
  model.par.c <- list("weibull" = 0.0, "median" = 0.3175, "hosking" = 0,
                      "blom" = 0.375, "cunnane" = 0.4, "gringorten" = 0.44,
                      "hazen" = 0.5)
  # model.types <- c("None", "normal", "lognormal", "exponential", "gumbel", "weibull",
  #                  "gev", "pearson3", "logpearson3", "uniform")
  # Normal: blom
  # Exponential + Gumbel: Gringorten
  # gev: cunnane
  # Pearson Type 3: Blom
  # model.par.c <- list("normal" = 0.375, "lognormal" = 0.375, "exponential" = 0.44,
  #                     "gumbel" = 0.44, "weibull" = 0.0, "gev" = 0.4, "pearson3" = 0.375,
  #                     "logpearson3" = 0.375, "uniform" = 0.0)
  Var1s <- sort(Var1, index.return = TRUE)
  n <- length(Var1)
  r <- 1:n
  #print((r-model.par.c[[model]]))
  Freq <- (r-model.par.c[[model]])/(n+1-2*model.par.c[[model]])
  results <- list(Var = Var1s$x, Freq = Freq, Treturn = 1.0/(1.0-Freq))
  return(results)
}
#' @title
#' probability_plot
#' @description
#' Function to calculate the measured and theoretical quantiles required in different probability
#' plots.
#' @param Var A numeric vector with the values of the variable to be checked
#' @param model A character string specifying the distribution function
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item Var: Measured variable
#' \item z: Theoretical variates
#' \item model: PDF model
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom lmom quape3
#' @importFrom lmom quagev
#' @importFrom stats qexp
#' @importFrom stats qunif
#' @importFrom stats qweibull
#' @export
probability_plot <- function(Var, model){
  Var1s <- sort(Var, index.return = TRUE)
  n <- length(Var)
  p <- (1:n)/(n+1)
  z <- NULL
  if(model == "normal"){
    z <- (p^.135-(1-p)^.135)/.1975
  }
  else if(model == "lognormal"){
    #Var1s$x <- log(Var1s$x)
    z <- (p^.135-(1-p)^.135)/.1975
  }
  else if(model == "exponential"){
    z <- qexp(p)
  }
  else if(model == "gumbel"){
    alpha <- sqrt(var(Var))*sqrt(6)/pi
    xi <- mean(Var) - 0.57772*alpha
    z <- as.vector(xi-alpha*log(-log(p)))
  }
  else if(model == "weibull"){
    z <- qweibull(p, shape = 1)
  }
  else if(model == "gev"){
    z <- quagev(p)
  }
  else if(model == "pearson3"){
    z <- quape3(p)
  }
  else if(model == "logpearson3"){
    #Var1s$x <- log(Var1s$x)
    z <- quape3(p)
  }
  else if(model == "uniform"){
    z <- qunif(p)
  }
  results <- list(Var = Var1s$x, z = z, model = model)
}
#' @title
#' calculate_diagram_moments
#' @description
#' Function to calculate the curves of the relationship between the Coefficient of Skewness and Kurtosis.
#' @return
#' A list with the following entries
#' \itemize{
#' \item N
#' \item LN3
#' \item P3
#' \item GEV
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
calculate_diagram_moments <- function(){
  Cs <- seq(-1.5, 4.0, 0.01)
  #three-parameter lognormal distribution (LN3):
  Ck <- 3 + 0.0256*Cs+1.7205*Cs^2+0.0417*Cs^3+0.0460*Cs^4-0.0048*Cs^5+0.0002*Cs^6
  LN3 <- cbind(Cs, Ck)
  #Pearson type III distribution (P3):
  Ck <- 3+1.5*Cs^2
  P3 <- cbind(Cs, Ck)
  # Generalized extreme-value distribution (GEV):
  Ck <- 2.6951 + 0.1858*Cs+1.7534*Cs^2+0.1107*Cs^3+0.0377*Cs^4+0.0036*Cs^5+
    0.0022*Cs^6 + 0.0007*Cs^7+0.00005*Cs^8
  GEV <- cbind(Cs, Ck)
  #
  results <- list(N = c(0.0,0.0), LN3 = LN3, P3 = P3, GEV = GEV)
  return(results)
}

