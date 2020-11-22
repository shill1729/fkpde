#' Solving tridiagonal linear systems fast
#'
#' @param alpha lower diagonal
#' @param beta diagonal
#' @param delta upper diagonal
#' @param b.vec the auxillary vectory in "Ax=b"
#'
#' @description {Thomas algorithm for solving tridiagonal systems fast without inverting any matrices.}
#' @details {The algorithm is fairly well known, i.e. on wikipedia and not complicated to implement. The tridiagonal system allows one to use back-substitution.}
#' @return vector
thomas_algorithm <- function(alpha, beta, delta, b.vec)
{
  n <- length(b.vec)
  cc <- matrix(0, nrow = n-1)
  dd <- matrix(0, nrow = n)
  x <- matrix(0, nrow = n)
  cc[1] <- delta[1]/beta[1]
  # Forward sweep for delta
  for(i in 2:(n-1))
  {
    cc[i] <- delta[i]/(beta[i]-alpha[i-1]*cc[i-1])
  }

  # Repeating the nearly same for b.vec
  dd[1] <- b.vec[1]/beta[1]
  for(i in 2:(n))
  {
    # alpha is given vector length of M-2, n = M-1. so i-1 ranges from 1 to n-1= 1 to M-2. looks good.
    dd[i] <- (b.vec[i]-alpha[i-1]*dd[i-1])/(beta[i]-alpha[i-1]*cc[i-1])
  }

  x[n] <- dd[n]
  for(i in (n-1):1)
  {
    x[i] <- dd[i]-cc[i]*x[i+1]
  }
  return(x)
}
