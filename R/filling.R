#' @title 
#' filling_average_station
#' @description 
#' Function to fill the missing data of a hydrologic record using the average of the station
#' @param Serie A numeric vector with the hydrologic time series
#' @return 
#' A numeric vector with the missing values filled with the station average
#' @author 
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family fill_data functions
#' @export
filling_average_station <- function(Serie){
  if(class(Serie) != "numeric"){
    stop("A numeric vector is required as input")
  }
  average_station <- mean(Serie, na.rm = TRUE)
  pos <- is.na(Serie)
  Serie1 <- Serie
  Serie1[pos] <- average_station
  return(Serie1)
}
#' @title 
#' filling_month_station
#' @description 
#' Function to fill the missing data of a hydrologic record using the monthly average
#' @param Serie A numeric vector with the hydrologic time series to be filled
#' @return 
#' A numeric vector with the corrected hydrologic time series
#' @author 
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family fill_data functions
#' @export
filling_month_station <- function(Serie){
  if(class(Serie) != "numeric"){
    stop("A numeric vector is required as input")
  }
  nyear <- floor(length(Serie) / 12)
  Serie1 <- pracma::Reshape(Serie, nyear, 12)
  average_month <- apply(Serie1, 1, mean, na.rm = TRUE)
  pos <-  is.na(Serie)

  return(Serie1)
}
#
filling_normal_ratio <- function(Serie1, Series){
  if(class(Serie1) != "numeric"){
    stop("A numeric vector is required as input")
  }
  #
  if(class(Series) != "matrix"){
    stop("A matrix is required as input")
  }
  #
  p0 <- mean(Serie1, na.rm = TRUE)
  pg <- apply(Series, 2, mean, na.rm = TRUE)
  weights <- p0/pg
  p0corr <- sum(weights*pg)
  return(p0corr)
}
#
filling_regression <- function(Serie1, Series){
  if(class(Serie1) != "numeric"){
    stop("A numeric vector is required as input")
  }
  #
  if(class(Series) != "matrix"){
    stop("A matrix is required as input")
  }
  #
  pos_valid <- !is.na(Serie1)
  Serie1_valid <- Serie1[pos_valid]
  Series_valid <- Series[pos_valid,]
  Series_valid1 <- cbind(Serie1_valid, Series_valid)
  Series.df <- as.data.frame(Series_valid1)
  #
  return(NULL)
}
