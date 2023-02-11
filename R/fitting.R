#' Fitting values for a tobit model
#'
#' This function calculates the different
#' fitted values from a probit model estimated
#' with AER package PyPos, EyPos, and Ey
#'
#' @param model The model to use
#' @return res: A matrix of the infile
#' @export
fitted.tobit = function(model){
  # input: an AER object
  # output: a 4xN matrix with initial values, PyPos EyPos and Ey

  # Check that model is an AER object
  if (! inherits(model, "survreg")) {
    stop("Error: object is not an AER object")
  }

  # Extract information
  x = model.matrix(model)
  b = model$coefficients
  s = model$scale
  val = model$y[, "time"]

  # Calculate the XB and XB/S
  xb = x%*% b
  xbs = xb/s

  # create PyPos
  PyPos = pnorm(xbs)

  # create EyPos (conditional expectations )
  EyPos =  xb + s * (dnorm(xbs)/pnorm(xbs))

  # create Ey (unconditional expectations)
  Ey = PyPos * EyPos

  res = cbind(val, PyPos, EyPos, Ey)
  colnames(res) = c("Value", "PyPos", "EyPos", "Ey")
  return(res)
}

