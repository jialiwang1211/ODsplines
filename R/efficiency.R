#' efficiency
#'
#' Compute the D-efficiency between two designs.
#'
#' @param t1 A vector of time points, usually the optimal design points.
#' @param t2 A vector of time points, usually the uniform design points.
#' @param eta smoothing parameter.
#' @param curvature square of prior curvature, function of t.
#'
#' @return A vector of optimal design points.
#'
#' @examples n<-4
#'drv2sq_fun<-curvature_logistic(theta1=1,theta2=-10,theta3=5)
#'OD_t<-ODsplines(n=n,curvature=drv2sq_fun)
#'t_u<-seq(0,1,length.out = n)
#'efficiency(OD_t,t_u,curvature=drv2sq_fun)
#'
#' @export
#' @import Deriv
#' @import expm
#' @import psych
efficiency<-function(t1,t2,eta=1,curvature){

obj_df<-function(t,eta,curvature,max){
    eta<-4/3*eta/max

    p<-length(t)

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

    obj<--log(det(t(X)%*%X)* det(t(Z_tilde)%*%Z_tilde+ eta*diag(p-2)- t(Z_tilde)%*%X%*%solve(t(X)%*%X)%*%t(X)%*%Z_tilde))
    df<-tr(solve(diag(p) + eta*delta%*%solve(Gnew)%*%t(delta)))
    return(list(obj=obj,df=df))
}

max<-max(curvature(seq(0,1,length.out = 100)))

a1<-exp(-obj_df(t1,eta,curvature,max)$obj)
df1<-obj_df(t1,eta,curvature,max)$df

a2<-exp(-obj_df(t2,eta,curvature,max)$obj)
df2<-obj_df(t2,eta,curvature,max)$df

return(a2^(1/df2)/a1^(1/df1))
}
