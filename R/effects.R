#' Extract the classes of the model independent variables
#' Class them into numeric and factor variables
varClass =  function(model){

  # extract the classes of each independent variables
  # except the intercept
  classes <- attributes(terms(model))[["dataClasses"]][-1L]

  # Check that model includes only numeric and facor variables
  # if not ask to recode some variables
  if (!all(classes %in% c("numeric", "factor", "logical"))) {
    stop("Error: The variables should be either numeric, logical, or factor")
  }

  fnames = NULL
  fn = unique(names(classes)[classes %in% c("factor")])

  if(! identical(fn, character(0))){
    fnames = paste0(unique(names(classes)[classes %in% c("factor")]), 1)}

  # identify factors versus numeric terms in `model`,
  #and examine only unique terms
  vars <- list(
    nnames = unique(names(classes)[classes %in% c("numeric")]),
    fnames = fnames
  )
  return(vars)
}

effect.cond.factor = function(X, b, s,  varname){
  # Reminder (EyPos =  xb + s * (dnorm(xbs)/pnorm(xbs)))

  X0 = X
  X0[, varname] <- 0
  xb0 = X0 %*% b
  xbs0 = xb0/s

  X1 = X
  X1[, varname] <- 1
  xb1 = X1 %*% b
  xbs1 = xb1/s

  E1Pos = xb1+ s *(dnorm(xbs1)/pnorm(xbs1))
  E0Pos = xb0+ s *(dnorm(xbs0)/pnorm(xbs0))

  res = cbind(mean(E1Pos-E0Pos), sd(E1Pos-E0Pos))
  return(res)
}

#' Calculate the average partial effects on conditional expected values
#'
#' This function calculates the different
#' fitted values from a probit model estimated
#' with AER package
#'
#' @param model The AER model to use
#'
#' @return
#' - call
#' - coefficients
#' - tableAPEcond
#'
#' @examples
#' library(wooldridge)
#' data(mroz)
#' model <- tobit(hours ~ nwifeinc + educ, data=mroz)
#' effect.cond(model)
#'
#' @export
effect.cond <- function(model){
  # Check that model is an AER object
  if (! inherits(model, "survreg")) {
    stop("Error: object is not an AER object")
  }

  # Extract information
  x = model.matrix(model)
  b = model$coefficients
  s = model$scale
  cname <- names(b)

  # Calculate the XB and XB/S
  xb = x%*% b
  xbs = xb/s

  # Extract the SE from AER model
  k_scale <- length(b)+1 # position of scale in thevar matrix
  vcov = model$var[-k_scale, -k_scale]
  colnames(vcov) = cname
  row.names(vcov) = cname
  stds = sqrt(diag(vcov))

  # Calculate the mean scaling factor to obtain
  # the marginal effects on conditionals for
  # continuous variables
  invM = dnorm(xbs)/pnorm(xbs)
  thet = 1 - invM * (xbs + invM)
  scale.cond = mean(thet)

  # Separate numerics from factors
  nnames = varClass(model)$nnames
  fnames = varClass(model)$fnames

  # Reporting of marginal effects on conditional
  # for continuous variables
  if (is.null(nnames)){
    table1 <- NULL
  } else{
  table1 <- matrix(rep(b[nnames], 4), ncol = 4)
  dimnames(table1) <- list(nnames, c("Value", "Std. Error", "z", "p"))
  table1[, 1] <- b[nnames] * scale.cond
  table1[, 2] <- stds[nnames] * scale.cond
  table1[, 3] <- table1[, 1]/stds[nnames]
  table1[, 4] <- 2*pnorm(-abs(table1[,3]))
  }

  # Reporting of marginal effects on conditional
  # for factors
  table2 <- NULL
  if (!is.null(fnames)){
  table2 <- matrix(rep(b[fnames], 4), ncol = 4)
  dimnames(table2) <- list(fnames, c("Value", "Std. Error", "z", "p"))
  for (i in fnames){
    table2[i, 1] <- effect.cond.factor(x,b,s, i)[,1]
    table2[i, 2] <- effect.cond.factor(x,b,s, i)[,2]
    table2[, 3]  <- table2[i, 1]/table2[i, 2]
    table2[, 4] <- 2*pnorm(-abs(table2[,3]))
  }
  }
  ## Combine tables
  if (is.null(fnames) & is.null(nnames)){
    table <- NULL
  } else {
  table <- rbind(table1, table2)
  }
  # returns object
  # return(list(table= table, tableN = table1, tableF=table2))
  res <- model[match(c('call', 'coefficients'),
                   names(model), nomatch=0)]
  res <- c(res, list(scaling = scale.cond, table=table,
                 tableN = table1, tableF = table2))
  class(res) <- 'tobitEffect'

  res
}

summary.tobitEffect <- function(x,
                                digits = max(options()$digits - 4, 3),
                                signif.stars=FALSE, ...) {

  # if(is.null(digits))
  #     digits <- options()$digits
  # cat("\nCall:\n")
  # dput(x$call)
  cat("\nAverage Partial Effects for Conditionals:\n")
  #    printCoefmat(x$table, digits = digits, signif.stars=signif.stars,
  #                 P.values=TRUE, has.Pvalue=TRUE)
  printCoefmat(x$table, digits = 3,
               has.Pvalue = TRUE, P.values = TRUE)
}

