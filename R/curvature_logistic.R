#' Square of curvature of logistic curve
#'
#' Compute the square of curvature values of 3-parameter logistic function.
#'
#' The 3-parameter logistic curve is given by \deqn{f(t)=\theta_1/(1+exp(\theta_2 t + \theta_3))} where \eqn{\theta_1} defines the asymptote, and \eqn{-\theta_3/\theta_2} defines the inflection point.
#'
#' @param theta1 parameter in logistic function that defines the asymptote.
#' @param theta2 parameter in logistic function.
#' @param theta3 parameter in logistic function.
#'
#' @return function of square of curvature values.
#'
#' @examples theta1<-1
#' theta2<--10
#' theta3<-5
#' t<-seq(0,1,by=0.01)
#' drv2sq_fun<-curvature_logistic(theta1,theta2,theta3)
#' plot(t,drv2sq_fun(t),type="l")
#'
#' @export
curvature_logistic<-function(theta1,theta2,theta3){

  curv<-function(t){
    x<-seq(0,1,by=0.01)

  f<-function(x){
  (-(theta1 * (exp(theta2 * x + theta3) * theta2 * theta2)/(1 +
  exp(theta2 * x + theta3))^2 - theta1 * (exp(theta2 * x +
  theta3) * theta2) * (2 * (exp(theta2 * x + theta3) * theta2 *
  (1 + exp(theta2 * x + theta3))))/((1 + exp(theta2 * x + theta3))^2)^2))^2
    }

  curvval<-rep(0,length(t))
  for (i in 1:(length(t)-1)){
    curvval[i]<-integrate(f, t[i], t[i+1])$value/(t[i+1]-t[i])
  }
  return(curvval)
  }

  return(curv)
}


