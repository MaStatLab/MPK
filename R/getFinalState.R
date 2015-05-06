#' Extract posterior draws from last iteration
#'
#'
#' @param object An \code{MPK} object
#' @param iter Integer indicating which iteration to extract. By default, the last iteration is extracted.
#' @return A list with the values of the MCMC at the \code{iter}-th iteration.
#' 
#' @examples
#' n = c(250, 250)
#' p = 4
#' 
#' Y1 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 3, ncol = p))
#' Y2 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 4, ncol = p))
#' Y = rbind(Y1, Y2)
#' C = c( rep(1,sum(n)), rep(2,sum(n)))
#' 
#' ans = mpk(Y, C)
#' state = getFinalState(ans)
#' ans2 = mpk(Y, C, state = state)


getFinalState <- function(object, iter = object$mcmc$nsave)
{
  state = list()  
  state$alpha = object$chain$alpha[iter]
  state$epsilon0 = object$chain$epsilon0[iter]
  state$epsilon = object$chain$epsilon[iter,]
  state$k_0 = object$chain$k_0[iter]
  state$m_1 = object$chain$m_1[iter,]
  state$mus = object$chain$mus[,,iter]
  state$mu_0 = object$chain$mu_0[,,iter]
  state$Omega_1 = object$chain$Omega_1[,,iter]
  state$Omegas = object$chain$Omegas[,,iter]
  state$Z = object$chain$Z[iter,]
  state$W = object$chain$w[,,iter]
  state$alphaHM = object$state$alphaHM
  state$epsilon0HM = object$state$epsilon0HM
  
  return(state)
}