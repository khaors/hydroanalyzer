#' @title
#' baseflow_graphical
#' @description
#' Function to calculate the baseflow from a discharge time series using the graphical approach
#' @param Q A numeric vector with the discharge time series
#' @return
#' A list with the following entries:
#' \itemize{
#' \item t: numeric vector with the time
#' \item Qbf: numeric vector with the calculated baseflow
#' }
#' @importFrom pracma interp1
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family graphical functions
#' @export
baseflow_graphical <- function(Q){
  if(class(Q) != 'numeric'){
    stop('A numeric vector is required as input')
  }
  dQs <- diff(Q)
  Qsmin <- vector('numeric', length=length(dQs))
  tmin <- vector('numeric', length=length(dQs))
  Qsmin[1:2] <- min(Q[1:2])
  pos <- 0
  for(i in 2:length(dQs)){
    p1 <- dQs[i-1]
    p2 <- dQs[i]
    if(p1 < 0.0 & p2 > 0.0){
      pos <- pos + 1
      Qsmin[pos] <- Q[i]
      tmin[pos] <- i
    }
  }
  #
  #print(tmin[1:pos])
  #print(Qsmin[1:pos])
  nt <- pos
  rt <- range(tmin[1:pos])
  bf <- interp1(tmin[1:nt],Qsmin[1:nt],rt[1]:rt[2])
  bf1 <- c(rep(bf[1],tmin[1]), bf, rep(bf[pos],length(Q)-rt[2]-1))
  #
  res <- list(t = 1:length(Q), Qbf = bf1)
  return(res)
}
