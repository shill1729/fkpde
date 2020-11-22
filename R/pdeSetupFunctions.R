#' Discretize function over space-time grid
#'
#' @param f function to evaluate on grid
#' @param x space variable
#' @param t time variable, can be null to discretize a function over just space.
#'
#' @description {Discretizes a function of space-time over a grid.}
#' @details {Simple double for-loop for space-time functions and one for-loop for functions of a single-variable}
#' @return matrix
gridFunction <- function(f, x, t = NULL)
{
  # If both x and t variable are present
  if(!is.null(t))
  {
    N <- length(t)
    M <- length(x)
    m <- matrix(data = 0, nrow = N, ncol = M)
    for(i in 1:N)
    {
      for(j in 1:M)
      {
        m[i, j] <- f(t[i], x[j])
      }
    }
  } else if(is.null(t)){
    n <- length(x)
    m <- matrix(data = 0, nrow = n)
    for(j in 1:n)
    {
      m[j] <- f(x[j])
    }
  }
  return(m)
}


#' The indicator function of an event
#'
#' @param bool a boolean for some event, statement, condition
#'
#' @description {A wrapper to \code{ifelse(bool, 1, 0)} to represent the indicator
#' of an event 'bool'.}
#' @details {Used like \code{indicator(x == 0)}, etc in probability theory.}
#' @return 1 if true, zero if false
#' @export indicator
indicator <- function(bool)
{
  ifelse(bool, 1, 0)
}

#' Initialize grids for the PDE solver
#'
#' @param dynamics list of functions defining the dynamics of the Ito process
#' @param problem list of functions defining the Feynman-Kac problem
#' @param region vector of maturity, lower and upper bounds of space
#' @param N time resolution
#' @param M space resolution
#' @description {Initializes grids used in the Feyman-Kac PDE solver.}
#' @details {The argument \code{dynamics} must be a named list of functions
#' \itemize{
#' \item \code{drift}, a function of \eqn{(t, x)} representing the local drift
#' \item \code{volat}, a function of \eqn{(t, x)} representing the local volatility
#' }
#' while the argument \code{problem} must be a named list of functions
#' \itemize{
#' \item \code{discount}, a function of \eqn{(t, x)} representing the discount function
#' \item \code{running_cost}, a function of \eqn{(t, x)} representing the running cost
#' \item \code{terminal_cost}, a function of \eqn{(x)} representing the terminal cost, as only a function of \eqn{x}.
#' }
#' and finally the argument \code{region} must be a vector with entries representing the maturity,
#' lower bound, upper bound, and time and space resolutions, in that order.}
#' @return list
#' @export setupPDE
setupPDE <- function(dynamics, problem, region, N, M)
{
  maturity <- region[1]
  A <- region[2]
  B <- region[3]


  grid.x <- seq(A, B, length.out = M+1)
  grid.t <- seq(0, maturity, length.out = N+1)
  h <- (B-A)/M
  k <- maturity/N
  # Grids, mu, volat.  1:N+1, 1:M+1; on paper 0:N, 0:M
  grids <- lapply(c(dynamics, problem[-3]), function(x) gridFunction(x, grid.x, maturity-grid.t))
  names(grids) <- c("m", "v", "r", "f")
  m <- grids$m
  v <- grids$v
  r <- grids$r
  f <- grids$f
  g <- gridFunction(problem[[3]], grid.x)

  alpha <- v^2/(2*h^2)-m/(2*h)
  beta <- -r-v^2/(h^2)
  delta <- v^2/(2*h^2)+m/(2*h)
  pde.setup <- list(grid.x = grid.x,
                    grid.t = grid.t,
                    alpha = alpha,
                    beta = beta,
                    delta = delta,
                    f = f,
                    g = g
  )
  return(pde.setup)
}



#' Compute conditional expectation via solving Feynman-Kac PDE
#'
#' @param dynamics list of functions defining the dynamics of the Ito process , see details
#' @param problem list of functions defining the Feynman-Kac problem, the conditional expectation
#' @param region vector, the space-time region to solve over and resolution of grids, see details
#' @param control list of time and space resolutions, variational boolean, and output option, and engine (r or c++)
#'
#' @description {Compute conditional expectations of functions of diffusions/Ito processes
#' via solving the Feynman-Kac PDE.}
#' @details {The argument \code{dynamics} must be a named list of functions
#' \itemize{
#' \item \code{drift}, a function of \eqn{(t, x)} representing the local drift
#' \item \code{volat}, a function of \eqn{(t, x)} representing the local volatility
#' }
#' while the arguments of \code{problem} must be a named list of functions
#' \itemize{
#' \item \code{discount}, a function of \eqn{(t, x)} representing the discount function
#' \item \code{running_cost}, a function of \eqn{(t, x)} representing the running cost
#' \item \code{terminal_cost}, a function of \eqn{(x)} representing the terminal cost, as only a function of \eqn{x}.
#' }
#' and finally the argument \code{region} must be a vector with entries representing the maturity,
#' lower bound, upper bound.}
#' @return list or numeric
#' @export solvePDE
solvePDE <- function(dynamics = NULL, problem = NULL, region = NULL, control = list(N = 100, M = 100, variational = TRUE, output = "greeks", engine = "c++"))
{

  N <- control$N
  M <- control$M
  variational <- control$variational
  output <- control$output
  engine <- control$engine
  if(is.null(dynamics))
  {
    stop("Need to pass function list defining model dynamics when initializing on solver call")
  }
  if(is.null(problem))
  {
    stop("Need to pass function list defining the conditional expectation when initializing on solver call")
  }
  if(is.null(region))
  {
    stop("Need to pass region/space-time resolution when initializing on solver call")
  }
  pde.setup <- setupPDE(dynamics, problem, region, N, M)
  h <- pde.setup$grid.x[2]-pde.setup$grid.x[1]
  k <- pde.setup$grid.t[2]-pde.setup$grid.t[1]
  g <- pde.setup$g
  f <- pde.setup$f
  u <- matrix(0, nrow = N+1, ncol = M+1)
  # IC
  u[1, ] <- g
  # BC
  u[,1] <- g[1]
  u[,M+1] <- g[M+1]

  # Create tridiagonal system
  for(i in 2:(N+1))
  {
    # auxiliary vector
    b.vec <- c(pde.setup$alpha[i, 2]*u[i, 1], rep(0, M-3), pde.setup$delta[i, M]*u[i, M+1])
    f.vec <- f[i-1, 2:M]
    b.vec <- b.vec + f.vec
    # Taking diagonals as vectors and using Thomas algorithm
    aa <- -k*pde.setup$alpha[i, 3:M]
    bb <- 1-k*pde.setup$beta[i, 2:M]
    cc <- -k*pde.setup$delta[i, 2:(M-1)]
    sol <- 0
    if(engine == "r")
    {
      sol <- thomas_algorithm(aa, bb, cc, k*b.vec+u[i-1, 2:M])
    } else if(engine == "c++")
    {
      sol <- thomasAlgorithm(aa, bb, cc, k*b.vec+u[i-1, 2:M])
    }
    if(variational){
      u[i, 2:M] <- pmax(sol, g[2:M])
    } else{
      u[i, 2:M] <- sol
    }

  }

  # Mixture volatility is undefined at zero, causes Naan
  # u <- u[stats::complete.cases(u), ]
  # N <- nrow(u)-1
  # M <- ncol(u)-1
  if(output == "grid")
  {
    return(list(space = pde.setup$grid.x,
                time = pde.setup$grid.t,
                solution = u
    )
    )
  } else if(output == "price")
  {
    return(u[N+1, M/2+1])
  } else if(output == "greeks")
  {
    i <- (M/2+1)
    fee <- u[N+1, i]
    delta <- (u[N+1, i+1]-u[N+1, i-1])/(2*h)
    gamma <- (u[N+1, i+1]-2*u[N+1, i]+u[N+1, i-1])/(h^2)
    # gamma <- (gamma-delta)/(spot^2)
    # delta <- delta/spot
    theta <- (u[N, i]-u[N+1, i])/k
    theta <- theta/360
    return(data.frame(fee = fee, delta = delta, gamma = gamma, theta = theta))
  }
}
