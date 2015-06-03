
library(MPK)
n = c(250, 250)
p = 4
 
Y1 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 3, ncol = p))
Y2 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 4, ncol = p))
Y = rbind(Y1, Y2)
C = c( rep(1,sum(n)), rep(2,sum(n)))

set.seed(1)
state = list(Z = sample(1:5, 1000, replace = TRUE))
mcmc = list(nburn = 1000, nsave = 1000, nskip = 1)
ans = mpk(Y, C, mcmc = mcmc, state = state)


cal = calibrate(ans)
cal.1 = calibrate(ans, reference.group = 1)
cal.2 = calibrate(ans, reference.group = 2)


par(mfrow=c(1,4))
plot(Y, col = C, xlim = range(Y[,1]), ylim = range(Y[,2]) )
plot(cal$Y_cal, col = C, xlim = range(Y[,1]), ylim = range(Y[,2]) )
plot(cal.1$Y_cal, col = C, xlim = range(Y[,1]), ylim = range(Y[,2]) )
plot(cal.2$Y_cal, col = C, xlim = range(Y[,1]), ylim = range(Y[,2]) )


plotDiff(ans, type = "weight")
plotDiff(ans, type = "shift")

plot.ts(ans$chain$epsilon0)

table(ans$chain$Z[1000,])

plot.ts(ans$chain$epsilon[,1])
plot.ts(ans$chain$epsilon[,4])

