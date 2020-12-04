# ----- functions for computing tail probabilities of truncated chi-squared distributions -----
# ----- full credit to Shuxiao Chen, the writer of the outference package  -----
#' A helper function for approximating normal tail probabilities
#'
#' For \eqn{Z ~ N(0, 1)}, we have the approximation
#'     \eqn{P(Z \ge z) \approx }\code{magicfun(z)*exp(-z^2/2)}.
#'
#' @keywords internal
#'
#' @param z, the number where the function is evaluated.
#'
#' @return This function returns the value of the function evaluated at \code{z}.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
magicfun = function(z){
  z2 <- z*z
  z3 <- z*z*z
  temp <- (z2 + 5.575192695 * z + 12.77436324) /
    (sqrt(2*pi) * z3 + 14.38718147*z2 + 31.53531977*z + 2*12.77436324)
  return(temp)
}

#' Make endpoints of intervals finite
#'
#' This function modifies a union of intervals with positive but possibly infinite endpoints
#'    into a union of intervals with positive and \emph{finite} endpoints, while ensuring
#'    the probability of a \eqn{N(0, 1)} falling into it numerically the same.
#'
#' @keywords internal
#'
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of intervals with \emph{positive} but possibly infinite endpoints.
#'
#' @return This function returns an "Intervals" object or a matrix depending on the input.
finiteE <- function(E) {
  ind.inf <- which(E == Inf)
  if (length(ind.inf) == 0) return(E)
  # we know there are some infinite entries
  E.max <- max(E[-ind.inf])
  E[which(E == Inf)] <- max(10000, E.max * 2)
  return(E)
}

#' Make endpoints of intervals positive
#'
#' This function modifies a union of intervals with possibly negative enpoints
#'     into a union of intervals with \emph{positive} endpoints, while ensuring
#'    the probability of a \eqn{N(0, 1)} falling into it numerically the same.
#'
#' @keywords internal
#'
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of intervals with \emph{positive} but possibly infinite endpoints.
#'
#' @return This function returns an "Intervals" object or a matrix depending on the input.
sortE <- function(E) {
  E.sorted <- lapply(1:nrow(E), function(i){
    temp <- as.numeric(E[i, ])
    if (temp[1] <= 0 & temp[2] <= 0) {
      return(sort(-temp))
    }
    if (temp[1] >= 0 & temp[2] >= 0) {
      return(sort(temp))
    }
    # we know temp[1] < 0, temp[2] > 0 OR temp[1] > 0, temp[2] < 0
    temp <- abs(temp)
    return(rbind(c(0, temp[1]), c(0, temp[2])))
  })
  E.sorted <- do.call(rbind, E.sorted)
  # in order to use the approximation, we translate Inf to a large number
  return(finiteE(E.sorted))
}

#' Comparison between two intervals
#'
#' This functions returns \code{TRUE} if and only if two intervals are the same.
#'
#' @keywords internal
#'
#' @param int1,int2 "Intervals" objects.
#'
#' @return This function returns the desired logical result.
isSameIntervals <- function(int1, int2) {

  # first make int1, int2 to the default order
  int1 <- intervals::reduce(int1)
  int2 <- intervals::reduce(int2)

  if (nrow(int1) != nrow(int2)) return(FALSE)

  # int1 and int2 has the same number of intervals

  if (sum(int1 != int2) > 0) return(FALSE)

  # int1 and int2 has the same elements
  return(TRUE)
}

#' Approximation of the ratio of two normal probabilities
#'
#' This function returns an approximation of \eqn{P(Z \in E1)/P(Z \in E2)}, where \eqn{Z ~ N(0, 1)}.
#'
#' @keywords internal
#'
#' @param E1,E2 "Intervals" objects or matrices where rows represents
#'     a union of intervals with \emph{positive and finite} endpoints.
#' @param scale scaling parameter.
#'
#' @return This function returns the value of the approximation.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
TNRatioApprox <- function(E1, E2, scale = NULL) {

  if (is.null(scale)) {
    temp <- (c(E1, E2))^2
    scale.grid <- stats::quantile(temp, probs = seq(0, 1, 0.2))

    for(scale in scale.grid) {
      temp <- TNRatioApprox(E1, E2, scale = scale)
      if (!is.na(temp)) {
        return(temp)
      }
      # if temp is NaN, proceed to the next loop
    }

    # if all scale.grid does not work, then return NaN
    return(NaN)
  }
  num1 <- magicfun(E1[, 1]) * exp(-(E1[, 1]^2 - scale)/2)
  num2 <- magicfun(E1[, 2]) * exp(-(E1[, 2]^2 - scale)/2)
  denom1 <- magicfun(E2[, 1]) * exp(-(E2[, 1]^2 - scale)/2)
  denom2 <- magicfun(E2[, 2]) * exp(-(E2[, 2]^2 - scale)/2)
  res <- sum(num1-num2)/sum(denom1-denom2)
  return(res)
}

#' Probability of a standard normal in a single interval
#'
#' This function returns \eqn{P(lo \le Z \le up)}, where \eqn{Z ~ N(0, 1)}.
#'
#' @keywords internal
#'
#' @param lo,up quantiles.
#'
#' @return This function returns the desired probability.
TNProbEachInt <- function(lo, up) {
  if (up == Inf) {
    return(stats::pnorm(lo, 0, 1, lower.tail = FALSE))
  }
  # we know up < Inf, want P(lo <= X <= up).
  # we want small - small (big - big will mask small numbers),

  try1 <- stats::pnorm(lo, 0, 1, lower.tail = FALSE) - stats::pnorm(up, 0, 1, lower.tail = FALSE)
  if (try1 != 0) return(try1)

  try2 <- stats::pnorm(up, 0, 1, lower.tail = TRUE) - stats::pnorm(lo, 0, 1, lower.tail = TRUE)
  return(try2)

}

#' Probability of a standard normal in a union of intervals
#'
#' This function returns \eqn{P(Z \in E)}, where \eqn{Z ~ N(0, 1)}.
#'
#' @keywords internal
#'
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of disjoint intervals.
#'
#' @return This function returns the desired probability.
TNProb <- function(E) {
  # sum cdf over each disjoint interval of E
  res <- sum(sapply(1:nrow(E), function(v) {
    return(TNProbEachInt(E[v, 1], E[v, 2]))
  }))
  return(res)
}

#' Survival function of truncated normal distribution
#'
#' This function returns the upper tail probability of a truncated normal distribution
#'     at quantile \code{q}.
#'
#' Let \eqn{X} be a normal random variable with \code{mean} and \code{sd}. Truncating
#'     \eqn{X} to the set \eqn{E} is equivalent to conditioning on \eqn{{X \in E}}. So this function
#'     returns \eqn{P(X \ge q | X \in E)}.
#'
#' @keywords internal
#'
#' @param q the quantile.
#' @param mean the mean parameter
#' @param sd the standard deviation
#' @param E the truncation set, an "Intervals" object or a matrix where rows represents
#'     a union of disjoint intervals.
#' @param approx should the approximation algorithm be used? Default is \code{FALSE},
#'     where the approximation is not used in the first place. But when the result is wacky,
#'     the approximation will be used.
#'
#' @return This function returns the value of the survival function evaluated at quantile \code{q}.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
TNSurv <- function(q, mean, sd, E, approx = FALSE) {
  # check if truncation is empty (i.e. truncated to the empty set)
  if (nrow(E) == 0) {
    stop("truncation set is empty")
  }

  # check if truncation is the whole real line
  if (isSameIntervals(E, intervals::Intervals(c(-Inf, Inf)))) {
    return(stats::pnorm(q, mean, sd, lower.tail = FALSE))
  }

  # E is not empty and is not the whole real line,
  # i.e. 0 < P(X in E) < 1

  # we want P(X > q | X in E) = P(X >= q AND X in E) / P(X in E)
  # {X >= q} = {Z >= (q-mean)/sd}
  # {X in E} = {Z in (E-mean)/sd}
  # Z ~ N(0, 1)
  q <- (q-mean)/sd
  E <- (E-mean)/sd
  mean <- 0
  sd <- 1
  q2 <- q*q
  region <- suppressWarnings(intervals::interval_intersection(E, intervals::Intervals(c(q, Inf))))
  # check if the result is 0 or 1
  if(nrow(region) == 0) return(0)
  if (isSameIntervals(E, region)) return(1)

  # transform region and E so that intervals have positive endpoints
  region <- sortE(region)
  E <- sortE(E)

  # we want P(Z in region) / P(Z in E)
  # try approximate calculation
  if (approx) {
    res <- TNRatioApprox(region, E, scale = q2)
    if (is.nan(res)) { # try approximation one more time
      res <- TNRatioApprox(region, E, scale = NULL)
    }
    return(max(0, min(1, res)))
  }

  # try exact calculation
  denom <- TNProb(E)
  num <- TNProb(region)

  if (denom < 1e-100 || num < 1e-100) {
    res <- TNRatioApprox(region, E, scale = q2)
    if (is.nan(res)) { # try approximation one more time
      res <- TNRatioApprox(region, E, scale = NULL)
    }
    return(max(0, min(1, res)))
  }

  # we know denom and num are both reasonably > 0

  res <- num / denom
  # force the result to lie in [0, 1]
  return(max(0, min(1, res)))
}

#' Approximation of the ratio of two chi-squared probabilities
#'
#' This function returns an approximation of \eqn{P(X \in E1)/P(X \in E2)}, where
#'     \eqn{X} is a central chi-squared random variable with \code{df} degrees of freedom.
#'
#' @keywords internal
#'
#' @param df degree of freedom of the chi-squared random variable.
#' @param E1,E2 "Intervals" objects or matrices where rows represents
#'     a union of intervals with \emph{positive and finite} endpoints.
#'
#' @return This function returns the value of the approximation.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
#' @references Canal, Luisa. "A normal approximation for the chi-square distribution."
#'     Computational statistics & data analysis 48.4 (2005): 803-808.
TChisqRatioApprox <- function(df, E1, E2) {

  # the transform that makes x into a N(0, 1) r.v. such that
  # P(X >= x) = P(Z >= Chisq2N(x)), X ~ chisq(df), Z ~ N(0, 1)
  # this function can take either scaler, vector or matrix
  Chisq2N <- function(x, df, tol = 1e-6) {

    if (is.numeric(x) && length(x) == 1) {
      if (x <= tol) { # x <= 0
        return(-Inf)
      }
      if (x == Inf) {
        return(Inf)
      }
      # we know x > 0 and x is finite
      x <- (x/df)^(1/6) - (1/2) * (x/df)^(1/3) + (1/3) * (x/df)^(1/2)
      mu <- 5/6 - 1/(9*df) - 7/(648*df^2) + 25/(2187*df^3)
      sig <- sqrt(1/(18*df) + 1/(162*df^2) - 37/(11664*df^3))
      return((x-mu)/sig)
    }

    if (is.vector(x)) {
      return(vapply(X = x, FUN = Chisq2N, FUN.VALUE = 0.1, df = df))
    }

    if (is.matrix(x)) {
      return(structure(vapply(X = x, FUN = Chisq2N, FUN.VALUE = 0.1, df = df), dim = dim(x)))
    }

    return(intervals::Intervals())

  }


  E1 <- Chisq2N(E1, df)
  E1 <- sortE(E1) # notice that Chisq2N can be negative
  E2 <- Chisq2N(E2, df)
  E2 <- sortE(E2)

  # now we want P(Z in E1) / P(Z in E2), Z ~ N(0, 1)
  return(TNRatioApprox(E1, E2))
}
