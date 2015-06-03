
setwd("~/Dropbox/Duke/Thesis/locally_tied_stick_breaking/examples/eqapol/30apr2015/")
library(MPK)

#####################################################
### same subject different stimulation conditions ###
### to see changes with stimulation               ###
#####################################################

data.1 = read.csv("~/Dropbox/Duke/Thesis/Flowcytometry/Equapol/EQAPOL_101_EP4_K69038GP_Unstim_02064.csv",  header = TRUE)
data.2 = read.csv("~/Dropbox/Duke/Thesis/Flowcytometry/Equapol/EQAPOL_101_EP4_K69038GP_CEF_02072.csv",  header = TRUE)
data.3 = read.csv("~/Dropbox/Duke/Thesis/Flowcytometry/Equapol/EQAPOL_101_EP4_K69038GP_CMVpp65_02069.csv",  header = TRUE)


remove.extremes <- function(data)
{
  for(i in 1:ncol(data))
  {
    idx.max = which( data[,i] == max(data[,i],na.rm = TRUE) )
    if(length(idx.max)>1)
    {
      data[idx.max, ] = NA
    }
    idx.min = which( data[,i] == min(data[,i],na.rm = TRUE) )
    if(length(idx.min)>1)
    {
      data[idx.min, ] = NA
    }  
  }    
  return(na.omit(data))
}


data.1 = remove.extremes(data.1)
data.2 = remove.extremes(data.2)
data.3 = remove.extremes(data.3)


data = rbind(data.1, data.2, data.3)


nobs = 5000
set.seed(1)
dim = seq(5,10)
Y_tot = scale( as.matrix( rbind(  data.1[sample(nrow(data.1),nobs),], 
                                  data.2[sample(nrow(data.2),nobs),],
                                  data.3[sample(nrow(data.3),nobs),] )  ) )[,dim]

C =  c(rep(1,nobs),rep(2,nobs), rep(3,nobs) )

idx.order = sample(nrow(Y_tot),nrow(Y_tot))
Y_tot = Y_tot[idx.order,]
C = C[idx.order]




mcmc = list(nburn = 5000, nsave = 500, nskip = 1, ndisplay = 100, seed = 2)

prior = list( K = 60 )

ans = mpk(Y_tot, C, prior = prior, mcmc = mcmc)



plotDiff(ans, type = "shift", dim = c(4,5))
plotDiff(ans, type = "weight", dim = c(4,5))

state = getFinalState(ans)

mcmc = list(nburn = 10, nsave = 500, nskip = 1, ndisplay = 100, seed = 2)

ans2 = mpk(Y_tot, C, prior = prior, mcmc = mcmc, state = state)

cal = calibrate(ans2)

par(mfrow = c(1,2))
plot(cal$Y_cal,  col = C+1)
plot(Y_tot,  col = C+1)
