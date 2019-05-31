#' optimal design for single population
#'
#' Derive the optimal design points of the adaptive smoothing splines. The prior knowledge is given  by the square of curvature values from a single population.
#'
#' @param n number of design points.
#' @param eta smoothing parameter.
#' @param curvature square of prior curvature, function of t.
#' @param interv a vector of length n-1, the mininal intervals between consecutive design points.
#' @param seed an integer.
#'
#' @return A vector of optimal design points.
#'
#' @examples n<-4
#'OD_t<-ODsplines(n=n,eta=1,
#'curvature=curvature_logistic(theta1=1,theta2=-10,theta3=5))
#'
#' @export
#' @import Deriv
#' @import expm
#' @import dfoptim
#' @import psych
ODsplines<-function(n,eta=1,curvature,interv=rep(0.001,n-1),seed=1){
  set.seed(seed)

  opt<-function(input){
    eta<-4/3*eta/max

    s<-input[-1]
    p<-length(s)
    t<-cumsum(s)

    a<-1.01
    b<-500

    h<-diff(t)

    lambda<-1/curvature(t)

    if(p==3){
      G<-(h[1:(p-2)]+h[2:(p-1)])/3
    }else{
      main_G<-(h[1:(p-2)]+h[2:(p-1)])/3
      upper_G<-lower_G<-h[2:(p-2)]/6
      G<-tridiag(upper=upper_G,lower=lower_G,main=main_G)
    }

    main_delta<-1/h[1:(p-2)]
    lower_delta<-1/h[2:(p-1)]
    middle_delta<--1/h[1:(p-2)]-1/h[2:(p-1)]
    delta<-bandtridiag(lower=lower_delta,middle=middle_delta,main=main_delta)

    if(p==3){
      Gstar<-(h[1:(p-2)]*lambda[1:(p-2)]+h[2:(p-1)]*lambda[2:(p-1)])/3
    }else{
      main_Gstar<-(h[1:(p-2)]*lambda[1:(p-2)]+h[2:(p-1)]*lambda[2:(p-1)])/3
      upper_Gstar=lower_Gstar<-h[2:(p-2)]*lambda[2:(p-2)]/6
      Gstar<-tridiag(upper=upper_Gstar,lower=lower_Gstar,main=main_Gstar)
    }

    X<-cbind(rep(1,p),t)
    Z<-delta%*%solve(t(delta)%*%delta)

    Gnew<-G%*%solve(Gstar)%*%G

    Z_tilde<-Z%*%sqrtm(Gnew)

    obj<--log(det(t(X)%*%X)* det(t(Z_tilde)%*%Z_tilde+ eta*diag(p-2)- t(Z_tilde)%*%X%*%solve(t(X)%*%X)%*%t(X)%*%Z_tilde))+
      exp((sum(s)-a)*b)
    return(obj)
  }

  t<-seq(0,1,length.out =n)
  p<-length(t)
  s<-t-c(0,t[-p])

  max<-max(curvature(seq(0,1,length.out = 500)))
  rs<-hjkb(c(eta,s), opt,
           lower=c(eta,c(0,interv)),
           upper=c(eta,rep(Inf,length(s))))

  OD_t<-cumsum(rs$par[-1])
  OD_t[length(t)]<-1


  return(OD_t)
}
