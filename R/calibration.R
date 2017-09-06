#' @title
#' sum_squared_residuals
#' @description
#' Function to calculate the sum of squared residuals between a vector of measured values and other
#' with the calculated values using a model.
#' @param measured A numeric vector with the measured values
#' @param calculated A numeric vector with the calculated values
#' @param ... Additional parameters
#' @return
#' A numeric value with the sum of squared residuals.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
sum_squared_residuals <- function(measured, calculated, ...){
  #return(sum((measured-calculated)^2))
  return(sqrt(mean((measured-calculated)^2)))
}
#' @title
#' mean_absolute_deviation
#' @description
#' Function to calculates the mean absolute deviation between a vector with measured values and
#' other vector with calculated values using a mathematical model.
#' @param measured A numeric vector with the measured values
#' @param calculated A numeric vector with the calculated values
#' @param ... Additional parameters
#' @return
#' A numeric value with the mean absolute deviation.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
mean_absolute_deviation <- function(measured, calculated, ...){
  return(mean(abs(calculated-measured)))
}
#' @title
#' max_absolute_deviation
#' @description
#' Function to calculates the maximum absolute deviation between a vector with measured values and
#' other vector with calculated values using a mathematical model.
#' @param measured A numeric vector with the measured values
#' @param calculated A numeric vector with the calculated values
#' @param ... Additional parameters
#' @return
#' A numeric value with the maximum absolute deviation.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
max_absolute_deviation <- function(measured, calculated, ...){
  return(max(abs(calculated-measured)))
}
#' @title
#' schultz_criteria
#' @description
#' Function to calculates the schultz criteria between a vector with measured values and
#' other vector with calculated values using a mathematical model. This criteria applies only
#' to dicharge data.
#' @param measured A numeric vector with the measured values
#' @param calculated A numeric vector with the calculated values
#' @param ... Additional parameters
#' @return
#' A numeric value with the schultz criteria.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
schultz_criteria <- function(measured, calculated, ...){
  num <- sum(abs(calculated-measured)*measured)
  den <- length(calculated)*(max(measured)^2)
  return(200*num/den)
}
#' @title
#' nash_sutcliffe
#' @description
#' Function to calculates the Nash-Sutcliffe criteria between a vector with measured values and
#' other vector with calculated values using a mathematical model. This criteria applies only
#' to dicharge data.
#' @param measured A numeric vector with the measured values
#' @param calculated A numeric vector with the calculated values
#' @param ... Additional parameters
#' @return
#' A numeric value with the  Nash-Sutcliffe criteria.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
nash_sutcliffe <- function(measured, calculated, ...){
  num <- sum((calculated-measured)^2)
  den <- sum((measured-mean(measured))^2)
  return(1.0-(num/den))
}
#' @title
#' log_nash_sutcliffe
#' @description
#' Function to calculates the logarithmic Nash-Sutcliffe criteria between a vector with measured values and
#' other vector with calculated values using a mathematical model. This criteria applies only
#' to dicharge data.
#' @param measured A numeric vector with the measured values
#' @param calculated A numeric vector with the calculated values
#' @param ... Additional parameters
#' @return
#' A numeric value with the logarithmic Nash-Sutcliffe criteria.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
log_nash_sutcliffe <- function(measured, calculated, ...){
  num <- sum((log(calculated)-log(measured))^2)
  den <- sum((log(measured)-mean(log(measured)))^2)
  return(1.0-(num/den))
}
#' @title
#' mass_balance_error
#' @description
#' Function to calculates the mass balance error between a vector with measured values and
#' other vector with calculated values using a mathematical model.
#' @param measured A numeric vector with the measured values
#' @param calculated A numeric vector with the calculated values
#' @param ... Additional parameters
#' @return
#' A numeric value with the mass balance error.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
mass_balance_error <- function(measured, calculated, ...){
  return(100*sum(calculated-measured)/sum(measured))
}

#' @title
#' objective.fn
#' @description
#' Calculates the objetive function
#' @param x A numeric vector with the parameters to be estimated
#' @param fn.model A character string with the name of the model used to evaluate the parameters
#' @param fit.fn A character string with the name of the function used to calculate the difference
#' between the calculated and measured values
#' @param params A list with the arguments of the fn.model
#' @param measured A numeric vector with the measured values
#' @return
#' A numeric value of the objective function
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
objective.fn <- function(x, fn.model, params, fit.fn, measured){
  args <- params
  args$par.model <- x
  calculated <- do.call(fn.model, args)
  #print(calculated)
  args1 <- list(measured = measured, calculated = calculated$Qt)
  res <- do.call(fit.fn, args1)
  return(res)
}
#' @title
#' calibrate
#' @description
#' Function to calibrate a water balance model using different optimization methods
#' @param measured A numeric vector with measured values
#' @param fn.model A character string with the name of the model function
#' @param obj.fn A character string with the name of the objective function
#' @param opt.method A character string with the name of the optimization method
#' @param args A list with a numeric vector with the precipitation values and a numeric
#' vector with the PEV values
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom stats optim
#' @importFrom GenSA GenSA
#' @importFrom GA ga
#' @export
calibrate <- function(measured, fn.model, obj.fn=c('rss', 'mnad', 'mxad'),
                      opt.method = c('lbfgs', 'ga', 'sa'), args){
  set.seed(12345)
  obj.fn1 <- ''
  par.res <- NULL
  res <- NULL
  if(obj.fn == 'rss'){
    obj.fn1 <- sum_squared_residuals
  }
  else if(obj.fn == 'mnad'){
    obj.fn1 <- mean_absolute_deviation
  }
  else if(obj.fn == 'mxad'){
    obj.fn1 <- max_absolute_deviation
  }
  #
  if(opt.method == 'lbfgs'){
    control.par <- list(maxit = 1000)
    #x <- c(.5, 400, .5, .1, 500, 50)
    x <- c(.5, 100, .5, .1, 100, 100)
    par.res <- optim(x, objective.fn, fn.model = fn.model, params = args, fit.fn = obj.fn1,
                     measured = measured, method = 'L-BFGS-B', lower = 1e-4,
                     upper = Inf, control = control.par)
    res <- list(par = par.res$par, value = par.res$value, opt.obj = par.res)
  }
  else if(opt.method == 'ga'){
    par.obj.fn <- function(x, fn.model, params, fit.fn, measured){
      -1*objective.fn(x, fn.model, params, fit.fn, measured)
    }
    lower <- rep(1e-4,6)
    upper <- c(1,1e3,1,.5,1e3,1e3)
    #control.par <- list(fnscale = -1)
    optimArgs.local <- list(method = "L-BFGS-B", potim = 0.2, pressel = 0.6)
    res.ga <- ga(type = "real-valued", fitness =  par.obj.fn, fn.model = fn.model,
                 params = args, fit.fn = obj.fn1, measured = measured,
                 min = lower, max = upper, popSize = 200, maxiter = 100, run = 50,
                 pcrossover = 0.8, pmutation = 0.1, parallel = FALSE,  optim = T,
                 optimArgs = optimArgs.local)
    res <- list( par = res.ga@solution, value = res.ga@fitnessValue, opt.obj = res.ga)
  }
  else if(opt.method == 'sa'){
    lower <- rep(1e-4,6)
    upper <- c(1,1e3,1,1,1e3,1e3)
    control.par <- list( verbose = TRUE, maxit = 50,
                         simple.function = T, max.time = 120)
    res.sa <- GenSA(lower = lower, upper = upper, fn = objective.fn, control = control.par,
                    fn.model = fn.model, params = args, fit.fn = obj.fn1,
                    measured = measured)
    res <- list(par = res.sa$par, value = res.sa$value, opt.obj = res.sa)
  }
  #
  return(res)
}
#' @title
#' uncertainty_quantification_boot
#' @description
#' Function to calculate the confidence intervals of the model parameters using a bootstrapping
#' approach.
#' @param measured A numeric vector with measured values
#' @param fn.model A character string with the name of the model function
#' @param obj.fn A character string with the name of the objective function
#' @param opt.method A character string with the name of the optimization method
#' @param args A list with a numeric vector with the precipitation values and a numeric
#' vector with the PEV values
#' @param level A numeric vector with the values of the quantiles associated with the
#' confidence interval to estimate
#' @param neval An integer specifying the number of model evaluation used to estimate the
#' confidence intervals
#' @param seed A random seed
#' @return
#' A list with the following entries:
#' \itemize{
#' \item
#' \item
#' }
#' @importFrom stats quantile
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
uncertainty_quantification_boot <- function(measured, fn.model,
                                            obj.fn=c('rss', 'mnad', 'mxad'),
                                            opt.method = c('lbfgs', 'ga', 'sa'), args,
                                            level = c(0.025, 0.975), neval = 100,
                                            seed = 12345){
  # Initial parameter estimation
  res.par <- calibrate(measured, fn.model, obj.fn = obj.fn, opt.method = opt.method,
                       args)
  npar <- length(res.par$par)
  # Calculate residuals
  fn.args <- args
  fn.args$par.model <- res.par$par
  fn.res <- do.call(fn.model, fn.args)
  calculated <- fn.res$Qt
  Q.residuals <- calculated - measured
  #
  set.seed(seed)
  #
  model.par <- matrix(0.0, nrow = neval, ncol = npar)
  for(ieval in 1:neval){
    current.residuals <- sample(Q.residuals, replace = TRUE)
    currentQ <- calculated + current.residuals
    current.res.par <- calibrate(currentQ, fn.model, obj.fn = obj.fn,
                                 opt.method = opt.method, args)
    model.par[ieval,] <- current.res.par$par
  }
  #
  par.ci <- matrix(0.0, nrow = npar, ncol = length(level))
  for(ipar in 1:npar){
    par.ci[ipar,] <- quantile(model.par[,ipar], probs = level)
  }
  #
  res <- list(model.par = model.par, par.ci = par.ci)
}
