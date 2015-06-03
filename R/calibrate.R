#' Calibration across data-sets
#'
#' Given draws from the posterior, this function calibrates data-sets.   
#'
#' @param x An \code{MPK} object
#' @param reference.group  
#' @return A list with three objects: 
#' \code{Ycal} is a matrix containing the calibrated data;
#' \code{calibration_distribution} is an array 
#' with the posterior distribution of the shift for each observation;
#' \code{calibration_median} is a matrix with the posterior median shift 
#' for each observation.
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
#' cal = calibrate(ans)
#' par(mfrow=c(1,2))
#' plot(Y, col = C)
#' plot(cal$Y_cal, col = C)
calibrate <- function(x, reference.group = NULL)
{
  if(class(x)!="MPK")
  {
    print("ERROR: Input's class should be MPK.")
    return(0)
  }
  C = x$data$C - 1
  Z = x$chain$Z - 1
  ref = ifelse(is.null(reference.group), -1, reference.group - 1 )
  output = calib(x$data$Y, 
                 C, 
                 Z, 
                 x$chain$mus, dim(x$chain$mus),
                 x$chain$mu_0, dim(x$chain$mu_0),
                 ref ) 
  colnames(output$Y_cal) = colnames(x$data$Y)
  output$reference.group = reference.group
  return(output)
  
}