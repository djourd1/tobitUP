#' Calculate the average partial effects on conditional expected values
#'
#' This function calculates the different
#' fitted values from a probit model estimated
#' with AER package
#'
#' @param model The model to use
#' @return res: A 4 columns matrix observed value, PyPos, EyPos, and Ey
#' @examples
#' library(wooldridge)
#' data(mroz)
#' model <- tobit(hours ~ nwifeinc + educ, data=mroz)
#' ape.tobit.cond(model)
#'
#' @export
ape.tobit.cond <- function(model){
  # input: an AER (survreg) object
  # output:

  # Check that model is an AER object
  if (! inherits(model, "survreg")) {
    stop("Error: object is not an AER object")
  }
  # Extract information
  x = model.matrix(model)
  b = model$coefficients

  cname <- names(model$coefficients)
  s = model$scale
  val = model$y[, "time"]
  kplus = length(b)+1

  # Calculate the XB and XB/S
  xb = x%*% b
  xbs = xb/s

  # thet = function(c) (1 - invMills(c) * (c + invMills(c)))
  invM = dnorm(xbs)/pnorm(xbs)
  thet = 1 - invM * (xbs + invM)
  scaleAPE.cond = mean(thet)

  # calculate APE assuming that all variables
  # ape.cond <- b * scaleAPE.cond
  table <- matrix(rep(b, 4), ncol = 4)
  dimnames(table) <- list(cname, c("Value", "Std. Error", "z", "p"))
  ape = b*scaleAPE.cond
  table[, 1] <- ape
  stds <- sqrt(diag(model$var)[-kplus])
  table[, 2] <- stds * scaleAPE.cond
  table[, 3] <- table[, 1]/stds
  table[, 4] <- 2*pnorm(-abs(table[,3]))
  x <- model[match(c('call', 'df', 'loglik', 'iter', 'na.action', 'idf',
                     'scale', 'coefficients'),
                   names(model), nomatch=0)]
  x <- c(x, list(scaleAPE.cond = scaleAPE.cond, tableAPEcond=table))
  class(x) <- 'tobitAPEcond'

  x
}

