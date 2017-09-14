
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
#
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


}
