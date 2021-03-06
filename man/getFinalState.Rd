% Generated by roxygen2 (4.0.2): do not edit by hand
\name{getFinalState}
\alias{getFinalState}
\title{Extract posterior draws from last iteration}
\usage{
getFinalState(object, iter = object$mcmc$nsave)
}
\arguments{
\item{object}{An \code{MPK} object}

\item{iter}{Integer indicating which iteration to extract. By default, the last iteration is extracted.}
}
\value{
A list with the values of the MCMC at the \code{iter}-th iteration.
}
\description{
Extract posterior draws from last iteration
}
\examples{
n = c(250, 250)
p = 4

Y1 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 3, ncol = p))
Y2 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 4, ncol = p))
Y = rbind(Y1, Y2)
C = c( rep(1,sum(n)), rep(2,sum(n)))

ans = mpk(Y, C)
state = getFinalState(ans)
ans2 = mpk(Y, C, state = state)
}

