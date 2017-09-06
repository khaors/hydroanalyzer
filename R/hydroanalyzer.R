#' HydroAnalyzer: Package to analyze hydrologic information
#'
#' @section calibration functions:
#'  The calibration functions includes functions to calibrate, calculate model fit, and objective
#'  functions.
#'
#'  The functions in this section are:
#'
#'  sum_squared_residuals, mean_average_deviation, nash_sutcliffe, log_nash_sutcliffe,
#'  mass_balance_error, calibrate, uncertainty_quantification
#'
#' @section digital_filter functions:
#'  The digital_filter functions include functions to calculate the baseflow from a discharge
#'  time series using different digital filters.
#'
#'  The functions in this section are:
#'
#'  nathan_filter, chapman_filter, eckhardt_filter
#'
#' @section graphical functions:
#'  The graphical functions include functions used to separate the baseflow from a discharge
#'  time series using simple graphical methods.
#'
#'  The function in this section is:
#'
#'  baseflow_graphical
#'
#' @section abcd functions:
#'  This section includes functions related to the abcd model. There are two version of this
#'  model implemented in HydroAnalyzer: monthly and yearly.
#'
#'  The functions in this section are:
#'
#'  abcd.year.model, abcd.month.model, abcd.mod.year.model, abcd.mod.month.model
#'
#' @section water budget funtions:
#' This sections includes functions related to the water budget models.
#'
#' The function in this section are:
#'
#' water_budget_direct
#'
#' @docType package
#' @name hydroanalyzer
NULL
