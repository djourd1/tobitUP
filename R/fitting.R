#' Fitting values for a tobit model
#'
#' This function calculates the different
#' fitted values from a probit model estimated
#' with AER package
#'
#' @param model The model to use
#' @return
#'
#' Returns a tobitExpected object with two values:
#'
#' + call: the tobit Call
#' + table: A 4 columns matrix observed value, PyPos, EyPos, and Ey
#'
#' @examples
#' library(wooldridge)
#' data(mroz)
#' model <- tobit(hours ~ nwifeinc + educ, data=mroz)
#' fitTobit(model)
#'
#' @export
fitTobit = function(model){
  # input: an AER object
  # output: a 4xN matrix with initial values, PyPos EyPos and Ey

  # Check that model is an AER object
  if (! inherits(model, "survreg")) {
    stop("Error: object is not an AER object")
  }

  # Extract information from the AER model
  x = model.matrix(model)
  b = model$coefficients
  s = model$scale
  val = model$y[, "time"]

  # Calculate the XB and XB/S
  xb = x%*% b
  xbs = xb/s

  # calculate PyPos
  PyPos = pnorm(xbs)

  # calculate EyPos (conditional expectations )
  EyPos =  xb + s * (dnorm(xbs)/pnorm(xbs))

  # Calculate Ey (unconditional expectations)
  Ey = PyPos * EyPos

  res = cbind(val, PyPos, EyPos, Ey)
  colnames(res) = c("Value", "PyPos", "EyPos", "Ey")

  res
}


#' Calculates the R2 of a tobit model
#'
#' This function calculates the pseudo R2
#' as the square of the correlation coefficient
#' between the observed values and the
#' unconditional fitted values E(y)
#'
#'
#' @param model The AER model
#' @return res: The R2 coefficient
#'
#' @examples
#' library(wooldridge)
#' data(mroz)
#' model <- tobit(hours ~ nwifeinc + educ, data=mroz)
#' r2Tobit(model)
#'
#' @export
r2Tobit = function(model){
  fits <- fitTobit(model)
  Ey = fits[, "Ey"]
  val = model$y[,"time"]
  r2<- cor(val, Ey)^2
  return(r2)
}



