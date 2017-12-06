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
  Var1s <- sort(Var1, index.return = TRUE)
  n <- length(Var1)
  r <- 1:n
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
  Var1 <- Var1s$x
  if(model == "normal"){
    z <- (p^.135-(1-p)^.135)/.1975
  }
  else if(model == "lognormal"){
    #Var1 <- log(Var1)
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
    #Var1 <- log(Var1)
    z <- quape3(p)
  }
  else if(model == "uniform"){
    z <- qunif(p)
  }
  results <- list(Var = Var1, z = z, model = model)
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
#' @export
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

calculate_parameters_moment <- function(x, pdf.model = "normal"){
  if(class(x) != numeric || class(x) != "matrix"){
    stop('ERROR: the required data is not numeric or matrix')
  }
  res <- NULL
  if(pdf.model == "normal"){
    res <- list(mean = mean(x), sd = sd(x))
  }
  else if(pdf.model == "lognormal"){
    mu <- mean(x)
    sigma <- sd(x)
    res <- list(mean = exp(mu+.5*sigma^2), sd = sqrt(exp(sigma^2)-1))
  }
  else if(pdf.model == "exponential"){
    res <- list(lambda = 1/mean(x))
  }
  else if(pdf.model == "gumbel"){
    beta <- (sqrt(6)/pi)*var(x)
    C <- -digamma(1)
    alpha <- mean(x)-beta*C
    res <- list(alpha = alpha, beta = beta)
  }
  else if(pdf.model == "weibull"){
    shape <- mean(x)
    scape <- 1
    res <- list(shape = shape, scale = scale)
  }
  else if(pdf.model == "gev"){

  }
  return(res)
}
#' @title
#' dgumbel
#' @description
#' Function to calculate the PDF of a Gumbel distribution
#' @param x A numeric vector with the values of the RV
#' @param mu The scale factor
#' @param s The shape factor
#' @return
#' A numeric vector with the density values of a Gumbel distribution.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
dgumbel <- function(x, mu, s){ # PDF
  exp((mu - x)/s - exp((mu - x)/s))/s
}
#' @title
#' pgumbel
#' @description
#' Function to calculate the CDF of a Gumbel distribution
#' @param q A numeric vector with quantiles of the RV
#' @param mu The scale factor
#' @param s The shape factor
#' @return
#' A numeric vector with the cumulative probability values of a Gumbel distribution.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
pgumbel <- function(q, mu, s){ # CDF
  exp(-exp(-((q - mu)/s)))
}
#' @title
#' qgumbel
#' @description
#' Function to calculate the quantiles of a Gumbel distribution
#' @param p A numeric vector with cumulative probability values
#' @param mu The scale factor
#' @param s The shape factor
#' @return
#' A numeric vector with the quantiles of a Gumbel distribution.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
qgumbel <- function(p, mu, s){ # quantile function
  mu-s*log(-log(p))
}
#' @title
#' estimate_parameters_mle
#' @description
#' Function to estimate the parameters of a given distribution using Maximum Likelihood
#' @param x A numeric vector with the RV values
#' @param pdf.model A character string with the name of the PDF to be used
#' @return
#' This function returns a list with the values of the estimated parameter for the PDF specified,
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom fitdistrplus fitdist
#' @importFrom stats sd
#' @export
estimate_parameters_mle <- function(x, pdf.model = "normal"){
  if(class(x) != "numeric"){
    stop('ERROR: the required data is not numeric')
  }
  res <- NULL
  if(pdf.model == "normal"){
    res <- list(mean = mean(x), sd = sd(x))
  }
  else if(pdf.model == "lognormal"){
    mu <- sum(log(x))/length(x)
    sigma <- sum((log(x)-mu)^2)/length(x)
    res <- list(mean = mu, sd = sqrt(sigma))
  }
  else if(pdf.model == "exponential"){
    lambda <- 1/mean(x)
    res <- list(lambda = lambda)
  }
  else if(pdf.model == "gumbel"){
    gumbel.fit <- fitdist(x, "gumbel", start=list(mu=5, s=5), method = "mle")
    res <- list(alpha = gumbel.fit$estimate[1], beta = gumbel.fit$estimate[2])
  }
  else if(pdf.model == "weibull"){
    weibull.fit <- fitdist(x, "weibull", start=list(shape=1, scale=1),
                           method = "mle")
    res <- list(shape = weibull.fit$estimate[1],
                scale = weibull.fit$estimate[2])
  }
  else if(pdf.model == "gev"){
    gev.fit <- fitdist(x, "gev", start = list(shape = 1, scale=1),
                       method = "mle")
    res <- list(shape = gev.fit$estimate[1],
                scale = gev.fit$estimate[2])

  }
  return(res)
}
