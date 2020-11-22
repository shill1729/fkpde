#' Animate the solution to a Feynman-Kac problem over time
#'
#' @param fk the solution grid returned from \code{solvePDE}
#' @param backwards boolean for terminal condition (TRUE) or initial condition (FALSE)
#' @param detail detail in \code{animate} from \code{gganimate}
#'
#' @description {Animates the time evolution of \eqn{u(t,x)} of the conditional expectation of a
#' function of a diffusion given its starting point.}
#' @return null
#' @export solutionAnimator
solutionAnimator <- function(fk, backwards = TRUE, detail = 10)
{
  N <- length(fk$time)-1
  M <- length(fk$space)-1
  if(backwards)
  {
    sol <- apply(fk$solution, 2, rev)
  } else
  {
    sol <- fk$solution
  }

  time_grid <- unlist(lapply(fk$time, function(x) rep(x, M+1)))
  dat <- data.frame(time = time_grid, price = rep(fk$space,N+1), value = as.numeric(t(sol)))
  # Animation
  ggplot2::aes_string()
  p <- ggplot2::ggplot(data = dat, ggplot2::aes_string("price", "value"))+
    ggplot2::geom_line()+
    gganimate::transition_time(dat$time)+
    gganimate::ease_aes('linear')+
    ggplot2::scale_x_continuous()+ggplot2::scale_y_continuous()+
    ggplot2::theme_minimal()
  ani <- gganimate::animate(p, end_pause = 10, detail = detail)
  print(ani)
  # Save time-stamped gif
  filename <- paste("fk ", Sys.time(), ".gif", sep = "")
  filename <- gsub(":", "-", filename)
  gganimate::anim_save(filename)

}
