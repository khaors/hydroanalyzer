#' @title 
#' von_neumann_test
#' @description 
#' Von Neumann test for homogenity of time series
#' @param Series A numeric vector with the hydrologic time series
#' @return 
#' A list with the following entries:
#' \itemize{
#' \item NN: A numeric value with the Von Neuman statistic
#' \item msg: A character string with the result of the test
#' }
#' @author 
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family homogenity functions
#' @export
von_neumann_test <- function(Series){
  nd <- length(Series)
  num <- sum((Series[1:nd-1]-Series[2:nd])**2)
  den <- sum((Series-mean(Series))**2)
  NN <- num/den
  msg <- ""
  #print("Von Neumann Test N= "%s+%NN)
  if(NN<2) {
    msg <- "Time series is non-homogeneous"
  }
  else {
    msg <- "Time series is homogeneous"
  }
  res <- list(NN = NN, msg = msg)
  return(res)
}
#' @title 
#' cumulative_deviation_test
#' @description 
#' Homogeneity test based on the adjusted partial sums or cumulative deviations from the mean following Buishand (1982).
#' @param Series A numeric vector with the hydrologic time series to be tested
#' @param prob The significance value to be used in the test
#' @return 
#' This function returns a list with the following entries:
#' \itemize{
#' \item Sk2: Rescaled adjusted partial sum 
#' \item Q: Change of level statistic
#' \item Qcrit: Critical value of Q
#' \item R: Range statistic
#' \item Rcrit: Critical value of R
#' \item msg: message
#' }
#' @importFrom stats sd
#' @export
#' @references 
#' Buishand, T. SOME METHODS FOR TESTING THE HOMOGENEITY OF RAINFALL RECORDS. Journal of 
#' Hydrology, 1982, 1, 11-27.
cumulative_deviation_test <- function(Series, prob){
  if(class(Series) != "numeric"){
    stop("ERROR: a numeric vector is required as input")
  }
  if(prob < 0.9){
    stop('ERROR: the specified confidence level is less than 0.9')
  }
  Sk1 <- cumsum(Series-mean(Series))
  Sk2 <- Sk1/sd(Series)
  Q <- max(abs(Sk2))
  R <- max(abs(Sk2))-min(abs(Sk1))
  ndat <- length(Series)
  pos_Q_n <- locate(table_Q_cumulative_deviation[,1],ndat)
  pos_R_n <- locate(table_R_cumulative_deviation[,1],ndat) 
  #
  if(abs(prob-0.9) < 1.e-10) {
    pos_Q_prob<-2
    pos_R_prob<-2
  }
  if(abs(prob-0.95) < 1.0e-10) {
    pos_Q_prob<-3
    pos_R_prob<-3
  }
  if(abs(prob-0.99) < 1.0e-10){
    pos_Q_prob<-4
    pos_R_prob<-4
  }
  Qtest <- table_Q_cumulative_deviation[pos_Q_n,pos_Q_prob]
  Rtest <- table_R_cumulative_deviation[pos_R_n,pos_R_prob]
  Q1<-Q/sqrt(ndat)
  R1<-R/sqrt(ndat)
  msg <- ""
  if(Q1<Qtest){
    msg <- "Time series is homogeneous"
  }
  else {
    msg <- "Time series is heterogeneous"
  }
  res <- list(Sk2 = Sk2, Q = Q, Qcrit = Qtest, R = R, Rcrit = Rtest, msg = msg)
  return(res)
}
#' @title 
#' locate
#' @description 
#' Function to locate a specific value inside a vector
#' @param vector_ A numeric vector
#' @param value a numeric value
#' @return 
#' This function returns the first position for which the value occurs in a vector
#' @author 
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family homogeneity functions 
locate <- function(vector_, value){
  ndat <- length(vector_)
  for(i in 2:ndat){
    if(vector_[i-1] <= value & vector_[i]>value) {
      return(i)  
    }
  }
  return(ndat)
}
