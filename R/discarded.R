#' Calculates the standard errors of the
#' partial effects on unconditional expected values
#' It requires a tobit model estimated with the  AER package
#'
#' @param model
#'
#' + model: The AER model to use
#' + R: the number of repetitions
#' + data: the data frame that was used for the model
#'
#' @return
#' - call: the call of the model
#' - table: A 4 x k table gathering the mean of
#' the average partial effets
#' - R: the number of repetitions
#'
#' @examples
#' library(wooldridge)
#' data(mroz)
#' model <- tobit(hours ~ nwifeinc + educ, data=mroz)
#' sem <- se.effect.uncond(model)
#' sem  #raw output
#' summary.uncond.eff(sem)  #output the formatted table
#'
#' @export
se.effect.uncond <- function(model, R, data, formu){
  effs <- NULL
  b = model$coefficients
  nr = nrow(data)
  print(nr)
  for (i in 1:R){
    new_sample = sample(1:nr, size = nr, replace = TRUE)
    df = data.frame(data[new_sample, ])
    new_model <- AER::tobit(formula = formu, data = df)
    effs = rbind(effs, effect.uncond(new_model)$vec  )
  }
  print(effs)
  table <- matrix(rep(b[-1], 4), ncol = 4)
  dimnames(table) <- list(names(b[-1]),
                          c("Value", "Std. Error", "z", "p")  )
  table[, 1] <- apply(effs, 2, mean)
  table[, 2] <- apply(effs, 2, sd)
  table[, 3] <- table[,1] / table[,2]
  table[, 4] <- 2*pnorm(-abs(table[,3]))

  # returns object
  res <- model[match(c('call', 'coefficients'),
                     names(model), nomatch=0)]
  res <- c(res, list(table=table, R = R))
  return(res)
}

#'
#' @export
summary.uncond.eff <- function(x,
                               digits = max(options()$digits - 4, 3),
                               signif.stars=TRUE, ...) {

  # if(is.null(digits))
  #     digits <- options()$digits
  # cat("\nCall:\n")
  # dput(x$call)
  cat("\nAverage Partial Effects for Unconditionals:\n")
  cat("(The standard errors calculating using bootstrap technique)\n\n")
  #    printCoefmat(x$table, digits = digits, signif.stars=signif.stars,
  #                 P.values=TRUE, has.Pvalue=TRUE)
  printCoefmat(x$table, digits = digits,
               has.Pvalue = TRUE,
               P.values = TRUE,
               signif.stars = signif.stars)
  cat("\nNo of repetitions: ")
  dput(x$R)

}
