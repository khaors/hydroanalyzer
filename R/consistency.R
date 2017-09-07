#' @title
#' double_mass_curve
#' @description
#' Function to check the consistency between two hydrological time series using the double mass curve approach.
#' @param Sref Numeric vector with the values of the reference time series
#' @param S1 Numeric vector with the values of the hydrological time series to be checked
#' @return
#' A list with the following entries:
#' \itemize{
#' \item S1ref: Cumulative values of the reference station
#' \item S1: Cumulative values of the testing stations
#' }
#' @export
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
double_mass_curve <- function(Sref, S1){
  ndat.ref <- length(Sref)
  ndat1 <- length(S1)
  #
  if(ndat.ref != ndat1){
    stop('ERROR: the number of data of the time series is not equal to the reference series')
  }
  #
  S1refs <- cumsum(Sref)
  S1s <- cumsum(S1)
  results <- list(S1ref = S1refs, S1 = S1s)
  return(results)
}
#' @title
#' bois_test
#' @description
#' This is a non-parametric test designed to check for the existence of a trend in the time series
#' of a single station or two stations using the cumulative residuals.
#' @param Serie1 A numeric vector with the values of the hydrologic time series to be checked
#' @param Serie2 A numeric vector with the values of the hydrologic time series to be checked
#' @param alpha Significance level of the test
#' @return
#' A list with the following entries:
#' \itemize{
#' \item residuals: Numeric vector with the residuals calculated for the testing station.
#' \item Ellipse.sup: Numeric vector with the values of the upper region of the confidence ellipse.
#' \item Ellipse.inf: Numeric vector with the values of the lower region of the confidence ellipse.
#' \item ylim: Numeric vector with the min and max value reached by the confidence ellipse.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom stats lm
#' @importFrom stats cor
#' @importFrom stats qnorm
#' @importFrom stats residuals
#' @importFrom stats var
#' @export
bois_test <- function(Serie1, Serie2, alpha){
  ndat1 <- length(Serie1)
  series.res <- NULL
  if(!missing(Serie2)){
    ndat2 <- length(Serie2)
    if(ndat1 != ndat2){
      stop('ERROR: the number of data in the time series is not equal')
    }
    serie.df <- data.frame(x = Serie1, y = Serie2)
    serie.lm <- lm(y ~ x , data = serie.df)
    sigma.y <- sqrt(var(serie.df$y))
    rho <- cor(serie.df$x, serie.df$y)
    sigma.res <- sigma.y*sqrt(1-rho^2)
    series.res <- residuals(serie.lm)
  }
  if(missing(Serie2)){
    series.res <- Serie1 - mean(Serie1)
    sigma.res <- sqrt(var(series.res))
  }
  E <- c(0, cumsum(series.res))
  jj <- 0:ndat1
  alpha1 <- 1.0-0.5*alpha
  z1 <- qnorm(alpha1)*sigma.res*sqrt((jj*(ndat1-jj)*(ndat1-1)/ndat1^2))
  z2 <- -1*qnorm(alpha1)*sigma.res*sqrt((jj*(ndat1-jj)*(ndat1-1)/ndat1^2))
  mn <- min(c(min(E),min(z2)))
  mx <- max(c(max(E),max(z1)))
  results <- list(residuals = E, ellipse.sup = z1, ellipse.inf = z2, ylim = c(mn,mx))
  return(results)
}
