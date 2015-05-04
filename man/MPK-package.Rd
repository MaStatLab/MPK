\name{MPK-package}
\alias{MPK-package}
\alias{MPK}
\docType{package}
\title{
Mixtures of Perturbed Kernels
}
\description{
Bayesian hierarchical model for comparison across mixtures distributions. 
         The model detects variations across the distributions in the mixture weights 
         and in the mixture locations.
}
\details{
\tabular{ll}{
Package: \tab MPK\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-04-29\cr
License: \tab BSD_3_clause + file LICENSE\cr
}

}
\author{
Jacopo Soriano
Maintainer: Jacopo Soriano <jacopo.soriano@gmail.com>
}
\references{
Literature or other references for background information
}

\keyword{ package }
\seealso{
 Optional links to other man pages, e.g. 
 \code{\link[<pkg>:<pkg>-package]{<pkg>}} 
}
\examples{
n = c(250, 250)
p = 4

Y1 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 3, ncol = p))
Y2 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 4, ncol = p))
Y = rbind(Y1, Y2)
C = c( rep(1,sum(n)), rep(2,sum(n)))

ans = mpk(Y, C)  
plotDiff(ans, type = "weight")
plotDiff(ans, type = "shift")
}
