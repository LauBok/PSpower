#' Calculate sample size needed to achieve a prespecified power
#'
#' @importFrom stats qnorm
#' @param tau the estimated treatment effect $E[Y(1) - Y(0)]$
#' @param alpha the significance level
#' @param beta the power to achieve
#' @param r the proportion of treated units
#' @param phi the overlap coefficients
#' @param E1,E0,S1,S0,R1,R0 the summary quantities
#' @param test whether one-sided or two-sided test is considered
#' @param estimand the estimand (ATE, ATT, ATC or ATO), or a customized tilting function
#' @returns an object with the calculated sample size
#' @examples
#' PSpower(1, 0.05, 0.956, 0.5, 0.99, -1.74, -2.74, 19.86, 20.12, 0.14, 0.14)
#' @export
PSpower <- function(tau, alpha, beta, r, phi, E1, E0, S1, S0, R1, R0,
                    test = 'two-sided', estimand = 'ATE') {
  params <- solve_parameters(r, phi, E1, E0, S1, S0, R1, R0)
  if (test == 'two-sided') {
    coef <- (qnorm(1 - alpha / 2) + qnorm(beta))^2
  }
  else if (test == 'one-sided') {
    coef <- (qnorm(1 - alpha) + qnorm(beta))^2
  }
  else {
    stop('test should be either one-sided or two-sided.')
  }
  if (typeof(estimand) == 'character') {
    if(!estimand %in% c('ATE', 'ATT', 'ATC', 'ATO'))
      stop('estimand should be one of ATE, ATT, ATC and ATO, or the tilting function.')
    if (estimand == 'ATE') {V <- params$V_ATE}
    else if (estimand == 'ATT') {V <- params$V_ATT}
    else if (estimand == 'ATC') {V <- params$V_ATC}
    else if (estimand == 'ATO') {V <- params$V_ATO}
  }
  else if (typeof(estimand) == 'closure') {
    V <- params$calc_V(estimand)
  }
  else {
    stop('estimand should be one of ATE, ATT, ATC and ATO, or the tilting function.')
  }
  return (structure(
    list(
      sample_size = coef * V / tau^2,
      tau = tau, alpha = alpha, beta = beta, test = test, estimand = estimand,
      summaries = c(r = r, phi = phi, E1 = E1, E0 = E0, S1 = S1, S0 = S0, R1 = R1, R0 = R0),
      params = params
    ),
    class = 'PSpower'
  ))
}

#' Prints PSpower object
#' @param x PSpower object
#' @param ... ignored
#' @export
print.PSpower <- function(x, ...) {
  cat('Estimated sample size:', x$sample_size, '\n')
}

#' Plots PSpower object
#' @param x PSpower object
#' @param power a range of powers to plot the power curve
#' @param ... ignored
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @examples
#' obj <- PSpower(1, 0.05, 0.956, 0.5, 0.99, -1.74, -2.74, 19.86, 20.12, 0.14, 0.14)
#' plot(obj)
#' @export
plot.PSpower <- function(x, power = seq(0.6, 0.99, length.out = 100), ...) {
  size <- sapply(power,
         function(p)
           PSpower(x$tau, x$alpha, p, x$summaries['r'], x$summaries['phi'],
                   x$summaries['E1'], x$summaries['E0'], x$summaries['S1'], x$summaries['S0'],
                   x$summaries['R1'], x$summaries['R0'], x$test, x$estimand)$sample_size
  )
  ggplot2::ggplot(data.frame(power = power, size = size), ggplot2::aes(size, power)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = data.frame(power = x$beta, size = x$sample_size))
}
