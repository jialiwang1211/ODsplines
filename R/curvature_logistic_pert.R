#' Square of curvature of logistic curve with perturbation
#'
#' Compute the square of curvature values of a mixture of 3-parameter logistic function and a Gaussian function.
#'
#' The curve is given by \deqn{f(t)=\theta_1/(1+exp(\theta_2 t + \theta_3)) - 0.02/\sigma * exp(-(t-\mu)^2/(2\sigma^2))}, where \eqn{\theta_1} defines the asymptote, and \eqn{-\theta_3/\theta_2} defines the inflection point.
#'
#' @param theta1 parameter in logistic function that defines the asymptote.
#' @param theta2 parameter in logistic function.
#' @param theta3 parameter in logistic function.
#' @param mu parameter in Gaussian function.
#' @param sigma parameter in Gaussian function.
#'
#' @return function of square of curvature values.
#'
#' @examples theta1<-1
#' theta2<--10
#' theta3<-5
#' mu<-0.7
#' sigma<-0.002
#' t<-seq(0,1,by=0.01)
#' drv2sq_fun<-curvature_logistic_pert(theta1,theta2,theta3,mu,sigma)
#' plot(t,drv2sq_fun(t),type="l")
#'
#' @export
curvature_logistic_pert<-function(theta1,theta2,theta3,mu,sigma){

curv<-function(t){
x<-seq(0,1,by=0.01)

f<-function(x){
(-(1 * theta1 * (exp(theta2 * x + theta3) * theta2 * theta2)/(1 +
exp(theta2 * x + theta3))^2 - 1 * theta1 * (exp(theta2 *
x + theta3) * theta2) * (2 * (exp(theta2 * x + theta3) *
theta2 * (1 + exp(theta2 * x + theta3))))/((1 + exp(theta2 *
x + theta3))^2)^2 + -0.02 * 1/sqrt(sigma * 2) * (exp(-(x -
mu)^2/(2 * sigma * 2)) * (2/(2 * sigma * 2)) - exp(-(x -
mu)^2/(2 * sigma * 2)) * (2 * (x - mu)/(2 * sigma * 2)) *
(2 * (x - mu)/(2 * sigma * 2)))))^2
}

curvval<-rep(0,length(t))
 for (i in 1:(length(t)-1)){
  curvval[i]<-integrate(f, t[i], t[i+1])$value/(t[i+1]-t[i])
  }
  return(curvval)
  }

  return(curv)
}


