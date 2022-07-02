# ----- general purpose helper functions -----

#' Takes the l2-norm of a vector.
#'
#' @keywords internal
#'
#' @param x the vector to be normed
#'
#' @return Returns the l2-norm of x.
norm_vec <- function(x) {
  sqrt(sum(x^2))
}

#' Checks if input is an integer between a and b
#'
#' @keywords internal
#'
#' @param x input to check 
#' @param a lower
#' @param b upper
#'
#' @return Returns TRUE if input is an integer between a and b, FALSE otherwise
is_integer_between_a_b <- function(x, a, b) {
  (x>= min(c(a, b))) && (x %% 1 == 0) && (x <= max(c(a, b)))
}

#' Checks if two clusterings are the same up to permutation
#'
#' @keywords internal
#'
#' @param cl1 the first clustering
#' @param cl2 the second clustering
#' @param K the number of clusters
#'
#' @return Returns TRUE if they are the same, and FALSE otherwise 
same_cl <- function(cl1, cl2, K) {
  tab <- table(cl1, cl2)
  sum(tab != 0) == K
}

#' Checks if Ck, Ck' in C(x'(phi))
#'
#' @keywords internal
#'
#' @param cl clustering of x
#' @param cl_phi clustering of x'(phi)
#' @param k1,k2 index of clusters involved in the test
#'
#' @return Returns TRUE if Ck, Ck' in C(x'(phi)), and FALSE otherwise
preserve_cl <- function(cl, cl_phi, k1, k2) {
  tab <- table(cl, cl_phi)

  k1_in <- (sum(tab[k1, ] != 0) == 1) & (sum(tab[, k1] != 0) == 1)
  k2_in <- (sum(tab[k2, ] != 0) == 1) & (sum(tab[, k2] != 0) == 1)

  k1_in & k2_in
}

#' Calculates log(gamma(q/2)) without overflowing float regardless of size of q
#'
#' @keywords internal
#'
#' @param q number of cols of dataframe
#'
#' @return Returns log(gamma(q/2))
calculate_log_gamma_q_over_2 <- function(q) {
  if (gamma(q / 2) != Inf) {
    return(log(gamma(q / 2)))
  }
  if (q %% 2 == 0) {
    return(calculate_log_gamma_q_even(q))
  }
  if (q %% 2 == 1) {
    return(calculate_log_gamma_q_odd(q))
  }
}

#' Calculates log(gamma(q/2)) if q > 343 (if gamma(q/2) overflows a float) and
#' q even.
#' Takes advantage of log rules (log(a*b) = log(a) + log(b)) and fact
#' that q will be positive integer.
#' Even function calculates log(gamma(q)) by summing values of log(x)
#' for(x in 1:((q/2)-1))
#' (because q/2 will be integer and gamma(integer) = (integer-1)!)
#'
#' @keywords internal
#'
#' @param q number of cols of dataframe (must be even)
#'
#' @return Returns log(gamma(q/2))
calculate_log_gamma_q_even <- function(q) {
  if (q %% 2 != 0) {
    stop("calculate_log_gamma_q_even can only take even values of q")
  }
  n <- (q / 2)
  # calculate log(n-1!)
  n_minus_one_fact_nums_to_sum <- c()
  for (x in 1:(n - 1)) {
    n_minus_one_factl_nums_to_sum <- c(n_minus_one_factl_nums_to_sum, log(x))
  }
  log_n_minus_one_fact <- sum(n_minus_one_fact_nums_to_sum)
  return(log_n_minus_one_fact)
}

#' Calculates log(gamma(q/2)) if q > 343 (if gamma(q/2) overflows a float) and
#' q odd
#' Takes advantage of log rules (log(a*b) = log(a) + log(b)),
#' (log(a/b) = log(a) - log(b)), and fact that q will be positive integer.
#'
#' @keywords internal
#'
#' @param q number of cols of dataframe (must be odd)
#'
#' @return Returns log(gamma(q/2))
calculate_log_gamma_q_odd <- function(q) {
  if (q %% 2 != 1) {
    stop("calculate_log_gamma_q_odd can only take odd values of q")
  }
  n <- ((q / 2) - .5)

  # calculate log((2n)! * sqrt(pi)), or equivalently log(2n!) + log(sqrt(pi))
  two_n_factorial_nums_to_sum <- c()
  for (x in 1:(2 * n)) {
    two_n_factorial_nums_to_sum <- c(two_n_factorial_nums_to_sum, log(x))
  }
  log_two_n_factorial <- sum(two_n_factorial_nums_to_sum)
  log_two_n_fact_sqrt_pi <- log_two_n_factorial + log(sqrt(pi))

  # calculate log((4^n) * n!), or equivalently log(4^n) + log(n!)
  log_four_to_power_of_n <- log(4) * n
  n_factorial_nums_to_sum <- c()
  for (x in 1:n) {
    n_factorial_nums_to_sum <- c(n_factorial_nums_to_sum, log(x))
  }
  log_n_factorial <- sum(n_factorial_nums_to_sum)
  log_four_to_n_n_factorial <- log_four_to_power_of_n + log_n_factorial

  # subtract first term (numerator) from second term (denominator)
  return(log_two_n_fact_sqrt_pi - log_four_to_n_n_factorial)
}
