#' Get differential score for each observation
#'
#'
#' @param object Object of class inheriting from "mpk".
#' @param type String indicating whether the score is associate to a difference in mixture weights
#' \code{type = "weight"}, or a difference in the mean parameters \code{type = "shift"}.
#' @return A vector with a score for each observation. 
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
#' score = getScore(ans, type = "shift")
#' hist(score)
getScore <- function(object, type = "weight")
{
  if(class(object)!="MPK")
  {
    print("ERROR: Input's class should be MPK.")
    return(0)
  }
  if(type == "weight"){
    state = apply(object$diff$logOdds, 2, median)    
  }else if(type == "shift"){
    state = apply(object$diff$shifts, 2, median)
  }else{
    print("ERROR: type is equal to 'weight' or 'shift'")
    return(0)
  }     
  return(state)
}