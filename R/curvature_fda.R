#' Square of curvature from "fd" object
#'
#' Estimate curvature using a functional data analysis approach, the curve is represented by basis functions, so that the curvature can be evaluated over the entire design region.
#'
#' @param basis a set of basis functions.
#' @param curvfd curvature function, "fd" class.
#'
#' @return function of square of curvature values.
#' @examples # Berkeley growth data for females
#' data(growth)
#' age<-growth$age
#' age<-(age-1)/17
#' heightmatf<-growth$hgtf
#'
#' #smoothing splines
#' norder<-6
#' nbasis<-length(age)+norder-2
#' heightbasis<-create.bspline.basis(c(0,2),nbasis, norder)
#' heightfdPar<-fdPar(heightbasis,3,1e-7)
#' heightfdSmoothf<-smooth.basis(age, heightmatf, heightfdPar)
#' heightfdf<-heightfdSmoothf$fd
#' accelfdUNf<-deriv.fd(heightfdf,2)
#' accelmeanfdUNf<-mean(accelfdUNf)
#'
#' #continuous registration
#' wbasisCR<-create.bspline.basis(c(0,2),15,6)
#' Wfd0CRf<-fd(matrix(0,15,dim(heightmatf)[2]),wbasisCR)
#' WfdParCRf<-fdPar(Wfd0CRf, 1, 1)
#' registerlistCRf = register.fd(accelmeanfdUNf, accelfdUNf , WfdParCRf)
#' accelfdCRf = registerlistCRf$regfd
#' accelmeanfdCRf = mean(accelfdCRf)
#'
#' drv2sq_fun<-curvature_fda(heightbasis,accelmeanfdCRf)
#'
#' @references
#' \insertRef{ramsay2005fda}{ODsplines}
#'
#' @export
#' @import fda
curvature_fda<-function(basis,curvfd){

  curv<-function(t){
    curvval<-rep(0,length(t))
    for (i in 1:(length(t)-1)){
      curvval[i]<- mean((predict(basis,seq(t[i],t[i+1],0.001))%*%curvfd$coefs)^2)

    }
    return(curvval)
  }

  return(curv)
}
