#' Fit Mixtures of Perturbed Kernels
#'
#' This function generates draws from the posterior of Mixtures of Perturbed Kernels. 
#'
#' @param Y Matrix of the data. Each row represents an observation.
#' @param C Vector of the group label of each observation. Labels are integers starting from 1. 
#' @param prior A list giving the prior information. If unspecified, a default prior is used. 
#' The list includes the following hyparameters: 
#' \code{K} Number of mixture components.
#' \code{merge_step} Introduce step to merge mixture components with small KL divergence. Default is \code{merge_step = TRUE}.
#' \code{merge_par} Parameter controlling merging radius. Default is \code{merge_par = 0.1}.
#' @param mcmc A list giving the MCMC parameters. If unspecified, default parameters are used.  
#' The list includes the following parameters: \code{nburn} indicates the number of burn-in scans,
#' \code{nsave} indicates the number of scans to be saved,
#' \code{nskip} indicates the thinning interval,
#' \code{ndisplay} indicates the number of scans to be displayed on screen. 
#' The total number of scans is equal to  \code{nburn + nsave*nskip}.
#' @param state Initial state of the chain.
#' @return A \code{MPK} object. 
#' @details
#' For \eqn{i = 1, \ldots, n_j} and \eqn{j = 1, \ldots, J}:
#' \deqn{y_{i,j} = \sum_{k=1}^{K}\pi_{j,k}N(y_{i,j} | \mu_{j,k}, \Sigma_k ),  }
#' where 
#' \deqn{(w_{j,1}, \ldots, w_{j,K}) | \alpha \sim Dirichlet(\alpha/K) }
#' \deqn{ \mu_{j,k} | \mu_{0,k} \epsilon_k, \Sigma_k \sim Normal(\mu_{0,k}, \epsilon_k \Sigma_k)}
#' \deqn{\epsilon_k | \epsilon_0  \sim Inv-Gamma( \tau_{\epsilon} + 1, \epsilon_0 \tau_{\epsilon}  ) }
#' \deqn{ \mu_{0,k} | \Sigma_k, k_0 \sim Normal(m_1, \Sigma_k / k_0)  }
#' \deqn{\Sigma_{k}^{-1}| \Psi_1 \sim Wishart(\Psi_1, \nu_1).}
#' In addition, there are the following hyperpriors:
#' \deqn{\alpha \sim Gamma(\tau_{\alpha, 1}, \tau_{\alpha, 2} ) }
#' \deqn{ \epsilon_0 \sim Beta(\tau_{\epsilon_0,1}, \tau_{\epsilon_0,2}  ) }
#' \deqn{ m_1 \sim Normal(m_2, S_2) } 
#' \deqn{k_0 \sim Gamma(\tau_{\gamma,1}/2, \tau_{\gamma,2}/2 )}
#' \deqn{ \Psi_1 \sim Inv-Wishart(\Psi_2, \nu_2)  }
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
mpk <- function(Y, C, prior = NULL, mcmc = NULL, state = NULL)
{
  Y = as.matrix(Y)  
  p = ncol(Y)
  
  J = length(unique(C))
  if( sum( sort(unique(C)) == 1:J )  != J )
  {
    print("ERROR: unique(C) should look like 1, 2, ...")
    return(0);
  }
  C = C - 1
  
  if(is.null(prior))
    prior = list()
  if(is.null(prior$K))
    prior$K = 10
  if(is.null(prior$m_2))
    prior$m_2 = colMeans(Y)
  if(is.null(prior$nu_2))
    prior$nu_2 = p+2
  if(is.null(prior$nu_1))
    prior$nu_1 = p+2
  if(is.null(prior$Psi_2))
    prior$Psi_2 = cov(Y)
  if(is.null(prior$S_2))
    prior$S_2 = 100*cov(Y)
  if(is.null(prior$tau_k0))
    prior$tau_k0 = c(4,4)
  if(is.null(prior$tau_alpha))
    prior$tau_alpha = c(1,1)
  if(is.null(prior$tau_epsilon))
    prior$tau_epsilon = p*J/2  
  if(is.null(prior$tau_epsilon0))
    prior$tau_epsilon0 = c(0.1,1) 
  if(is.null(prior$merge_step))
    prior$merge_step = TRUE
  if(is.null(prior$merge_par))
    prior$merge_par = 0.1
    
      
  if(is.null(mcmc))
    mcmc = list()
  if(is.null(mcmc$nburn))
    mcmc$nburn = 5000
  if(is.null(mcmc$nsave))
    mcmc$nsave = 1000
  if(is.null(mcmc$nskip))
    mcmc$nskip = 1
  if(is.null(mcmc$ndisplay))
    mcmc$ndisplay = 100
  if(is.null(mcmc$seed))
    mcmc$seed = 42
    

  # set initial random seed for R here, it will propagate down to Rcpp
  # but Armadillo's RNG will get set in MCMC class in gibbs.h
  # NOTE: Rstudio seems to touch Armadillo's underlying RNG, so
  # differing results may be seen when run from Rstudio, see comments
  # from Dirk Eddelbuettel here:
  #  https://github.com/RcppCore/RcppArmadillo/issues/11
  set.seed(mcmc$seed) 

  
  if(is.null(state))
    state = list()  
  
  if(is.null(state$Z))
    state$Z = kmeans(Y, 2*prior$K/3, iter.max = 100, algorithm="Lloyd")$cluster - 1
  
  if(is.null(state$epsilon))
    state$epsilon = 1/rgamma(prior$K, shape = prior$tau_epsilon, 
                             rate = prior$tau_epsilon*prior$tau_epsilon0[1]/sum(prior$tau_epsilon0))
  
  if(is.null(state$epsilon_0))
    state$epsilon_0 = prior$tau_epsilon0[1]/sum(prior$tau_epsilon0)
  
  if(is.null(state$Omega_1))
    state$Omega_1 = solve(cov(Y))
  
  if(is.null(state$Omegas))
  {    
    N = tabulate(as.numeric(state$Z), prior$K)
    state$Omegas = matrix(NA, nrow = p , ncol = p*prior$K)
    for(k in 1:prior$K)
    {
      if(N[k] > p + 1)
        state$Omegas[,seq(p*(k-1)+1,p*(k-1)+p)] = solve(cov(Y[state$Z==k,]))
      else
        state$Omegas[,seq(p*(k-1)+1,p*(k-1)+p)] = state$Omega_1   
    }
  }
  
  if(is.null(state$mu_0))
  {
    N = tabulate(as.numeric(state$Z), prior$K)
    state$mu_0 = matrix(NA,nrow=p, ncol=prior$K)
    for(k in 1:prior$K)
    {
      if(N[k] > 0)
        state$mu_0[,k] = colMeans(Y[state$Z==k,]) 
      else
        state$mu_0[,k] = colMeans(Y) + rnorm(p)  # do better if possible
    }
  }  
  
  if(is.null(state$mus))
  {
    state$mus = matrix(NA,nrow=J, ncol=p*prior$K)
    for(k in 1:prior$K)
    {
      for(j in 1:J)
        state$mus[j,seq(p*(k-1)+1,p*(k-1)+p)] = state$mu_0[,k] + rnorm(p) # do better if possible 
    }

  }  
  
  if(is.null(state$alpha))
    state$alpha = prior$tau_alpha[1]/sum(prior$tau_alpha)
  
  if(is.null(state$k_0))
    state$k_0 = prior$tau_k0[1]/sum(prior$tau_k0)

  if(is.null(state$w))
    state$w = matrix(1/prior$K, nrow = J, ncol = prior$K)
  
  if(is.null(state$m_1))
     state$m_1 = colMeans(Y)  

  #####################################
  
  ans = fit(Y, C, prior, mcmc, state)
  colnames(ans$data$Y) = colnames(Y)
  ans$data$C = ans$data$C + 1
  ans$chain$Z = ans$chain$Z + 1 
  class(ans) = "MPK"
  return(ans)
  
}
