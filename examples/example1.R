
library(MPK)
n = c(150, 150)
p = 4
 
Y1 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 3, ncol = p))
Y2 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 4, ncol = p))
Y = rbind(Y1, Y2)
C = c( rep(1,sum(n)), rep(2,sum(n)))
 
ans = mpk(Y, C)

plotDiff(ans, type = "weight")
plotDiff(ans, type = "shift")



