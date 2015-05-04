#' Plot differences across mixture models
#'
#'
#' @param object Object of class inheriting from "mpk".
#' @param dim Vector indicating along which dimension the data are projected. 
#' Default is \code{dim = c(1,2)}.
#' @param type String indicating whether the data-points are colored according to difference in mixture weights
#' \code{type = "weight"}, or difference in the mean parameters \code{type = "shift"}.
#' @param group Vector indicating which groups are plotted. By default all groups are plotted.
#' @param xlim the x limits (x1, x2) of the plot.
#' @param xlim the y limits of the plot.
#' @return A plot. 
#' @examples
#' 
#' n = c(250, 250)
#' p = 4
#' 
#' Y1 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 3, ncol = p))
#' Y2 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 4, ncol = p))
#' Y = rbind(Y1, Y2)
#' C = c( rep(1,sum(n)), rep(2,sum(n)))
#' 
#' ans = mpk(Y, C)
#' plotDiff(ans, type = "weight")
#' plotDiff(ans, type = "shift")
plotDiff <- function( object, type = "weight", 
                      dim = c(1,2), group =  1:length(unique(object$data$C)),
                      xlim = range(object$data$Y[,dim[1]]), 
                      ylim = range(object$data$Y[,dim[2]]), ...  )
{
  if(class(object)!="MPK")
  {
    print("ERROR: Input's class should be MPK.")
    return(0)
  }
  
  colors = rainbow(150)[1:100]
  if(type == "weight"){
    state = apply(object$diff$logOdds, 2, median)    
  }else if(type == "shift"){
    state = apply(object$diff$shifts, 2, median)
  }else{
    print("ERROR: type is equal to 'weight' or 'shift'")
    return(0)
  }     
  
  
  size = ceiling(sqrt(length(group)))
  if( size * (size -1) == length(group) ){
    mat = matrix( seq(1,  size*(size-1)) , ncol=size, byrow=TRUE  )
  }else{
    mat = matrix( seq(1,  size*size) , ncol=size, byrow=TRUE  ) 
  }
  layout(mat = mat)
  
  if(is.null(colnames(object$data$Y)[dim[1]])){
    xlab = ""
  }else
  {
    xlab = colnames(object$data$Y)[dim[1]]
  }
  
  if(is.null(colnames(object$data$Y)[dim[2]]))
    ylab = ""
  else
    ylab = colnames(object$data$Y)[dim[2]]
  
  for(j in group)
  {    
    temp = state[object$data$C==j]
    if(type == "shift"){      
      effect = round(temp/max(max(state),1)*99) + 1
      order = sort.int(temp, index.return = TRUE)$ix
    }else{
      effect = round((temp - min(-1,min(state))) / (max(max(state),1) - min(-1,min(state))) *99) + 1
      order = sort.int(abs(temp), index.return = TRUE)$ix
    }

    plot(object$data$Y[object$data$C==j,dim][order,],
         col = colors[effect[order]], 
         pch = 16, xlim = xlim*c(1,1.1), ylim = ylim,
         xlab = xlab,
         ylab = ylab)
    par(new = TRUE)
    plot(rep(10,100), (1:100)/100, col = colors[1:100], bty='l', axes = FALSE,
         pch = 15, xlab = "", ylab = "", xaxt='n', cex = 2, 
         ylim = c(0,1), xlim = c(0,10))
    if(type == "shift")
      labels = signif(seq( 0, max(max(state),1),length.out=6), digits=2)
    else
      labels = signif(seq( min(-1,min(state)), max(max(state),1),length.out=6), digits=2)
    axis(side = 4, at=seq(0,1,by = .2), labels=labels, las = 2)
    par(new = FALSE)
  }  
  
  layout(1)  
  
}