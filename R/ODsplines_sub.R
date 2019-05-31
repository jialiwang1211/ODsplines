#' Optimal design points for two populations
#'
#' Derive the optimal design points of the adaptive smoothing splines with different prior curvatures from two subpopulation.
#'
#' @param n number of design points.
#' @param eta smoothing parameter.
#' @param curvature1 square of prior curvature of population 1, function of t.
#' @param curvature2 square of prior curvature of population 2, function of t.
#' @param interv a vector of length n-1, the mininal intervals between consecutive design points.
#' @param seed an integer.
#'
#' @examples n<-7
#' curvature1<-curvature_logistic(theta1=1,theta2=-12,theta3=5)
#' curvature2<-curvature_logistic(theta1=1,theta2=-9,theta3=6)
#' OD_t<-ODsplines_sub(n=n,curvature1=curvature1,curvature2=curvature2)
#'
#' @export
#' @import Deriv
#' @import expm
#' @import dfoptim
#' @import psych
#' @import Matrix
ODsplines_sub<-function(n,eta=1,curvature1,curvature2,interv=rep(0.001,n-1),seed=1){
set.seed(seed)

opt<-function(input){
eta<-4/3*eta/max

s<-input[-1]
p<-length(s)
t<-cumsum(s)

a<-1.01
b<-500

h<-diff(t)

lambda1<-1/curvature1(t)
lambda2<-1/curvature2(t)

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
  Gstar1<-(h[1:(p-2)]*lambda1[1:(p-2)]+h[2:(p-1)]*lambda1[2:(p-1)])/3
  Gstar2<-(h[1:(p-2)]*lambda2[1:(p-2)]+h[2:(p-1)]*lambda2[2:(p-1)])/3
}else{
  main_Gstar1<-(h[1:(p-2)]*lambda1[1:(p-2)]+h[2:(p-1)]*lambda1[2:(p-1)])/3
  upper_Gstar1=lower_Gstar1<-h[2:(p-2)]*lambda1[2:(p-2)]/6
  main_Gstar2<-(h[1:(p-2)]*lambda2[1:(p-2)]+h[2:(p-1)]*lambda2[2:(p-1)])/3
  upper_Gstar2=lower_Gstar2<-h[2:(p-2)]*lambda2[2:(p-2)]/6

  Gstar1<-tridiag(upper=upper_Gstar1,lower=lower_Gstar1,main=main_Gstar1)
  Gstar2<-tridiag(upper=upper_Gstar2,lower=lower_Gstar2,main=main_Gstar2)
}

X<-cbind(rep(1,p),t)
Z<-delta%*%solve(t(delta)%*%delta)

Gnew1<-G%*%solve(Gstar1)%*%G
Gnew2<-G%*%solve(Gstar2)%*%G

Z_tilde1<-Z%*%sqrtm(Gnew1)
Z_tilde2<-Z%*%sqrtm(Gnew2)

bigX<-rbind(X,X)
bigZ<-as.matrix(bdiag(Z_tilde1,Z_tilde2))

obj<--log(det(t(bigX)%*%bigX)* det(t(bigZ)%*%bigZ+ eta*diag(2*(p-2))- t(bigZ)%*%bigX%*%solve(t(bigX)%*%bigX)%*%t(bigX)%*%bigZ))+
  exp((sum(s)-a)*b)

return(obj)
}

t<-seq(0,1,length.out =n)
p<-length(t)
s<-t-c(0,t[-p])

max<-max(c(curvature1(seq(0,1,length.out = 500)),curvature2(seq(0,1,length.out = 500))))
rs<-hjkb(c(eta,s), opt,
         lower=c(eta,c(0,interv)),
         upper=c(eta,rep(Inf,length(s))))

OD_t<-cumsum(rs$par[-1])
OD_t[length(t)]<-1
return(OD_t)
}

