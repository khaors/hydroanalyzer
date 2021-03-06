% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration.R
\name{schultz_criteria}
\alias{schultz_criteria}
\title{schultz_criteria}
\usage{
schultz_criteria(measured, calculated, ...)
}
\arguments{
\item{measured}{A numeric vector with the measured values}

\item{calculated}{A numeric vector with the calculated values}

\item{...}{Additional parameters}
}
\value{
A numeric value with the schultz criteria.
}
\description{
Function to calculates the schultz criteria between a vector with measured values and
other vector with calculated values using a mathematical model. This criteria applies only
to dicharge data.
}
\seealso{
Other calibration functions: \code{\link{log_nash_sutcliffe}},
  \code{\link{mass_balance_error}},
  \code{\link{max_absolute_deviation}},
  \code{\link{mean_absolute_deviation}},
  \code{\link{nash_sutcliffe}}, \code{\link{objective.fn}},
  \code{\link{sum_squared_residuals}},
  \code{\link{uncertainty_quantification_boot}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
\concept{calibration functions}
