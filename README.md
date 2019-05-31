# ODsplines

Derive the D-optimal design points for temporal data using adaptive smoothing splines, where the prior knowledge is the curvature values. \
The function 'ODsplines' performs optimal designs for a single population. The argument 'curvature' can be provided from a parametric function, for example the curvature of a logistic curve using the function 'curvature_logistic', or a mixture of logistic curve and a Gaussian function using the function 'curvature_logistic_pert'. \
The function 'ODsplines_sub' performs optimal designs for two populations with different curvature values, which should be provided from the arguments 'curvature1' and 'curvature2'. \
Application examples to the Berkeley growth data are illustrated below for 1) females, and 2) both males and females.  


### Installation: 
library(devtools)\
install_github("jialiwang1211/ODsplines")\
library(ODsplines)

### Berkeley growth data examples
library(fda)\
data(growth)\
age<-growth$age\
age<-(age-1)/17\
heightmatf<-growth$hgtf

### #1. females
#### #smoothing splines
norder<-6\
nbasis<-length(age)+norder-2\
heightbasis<-create.bspline.basis(c(0,2),nbasis,norder)\
heightfdPar<-fdPar(heightbasis,3,1e-7)\
heightfdSmoothf<-smooth.basis(age,heightmatf,heightfdPar)\
heightfdf<-heightfdSmoothf$fd\
accelfdUNf<-deriv.fd(heightfdf,2)\
accelmeanfdUNf<-mean(accelfdUNf)

#### #continuous registration
wbasisCR<-create.bspline.basis(c(0,2),15,6)\
Wfd0CRf<-fd(matrix(0,15,dim(heightmatf)[2]),wbasisCR)\
WfdParCRf<-fdPar(Wfd0CRf,1,1)\
registerlistCRf<-register.fd(accelmeanfdUNf,accelfdUNf,WfdParCRf)\
accelfdCRf<-registerlistCRf$regfd\
accelmeanfdCRf<-mean(accelfdCRf)

drv2sq_funf<-curvature_fda(heightbasis,accelmeanfdCRf)

#### #OD for females using function 'ODsplines'
n<-30\
OD_f<-ODsplines(n=n,curvature=drv2sq_funf)

### #2. both females and males
#### #repeat for males
heightmatm<-growth$hgtm\
heightfdSmoothm<-smooth.basis(age,heightmatm,heightfdPar)\
heightfdm<-heightfdSmoothm$fd\
accelfdUNm<-deriv.fd(heightfdm,2)\
accelmeanfdUNm<-mean(accelfdUNm)\
Wfd0CRm<-fd(matrix(0,15,dim(heightmatm)[2]),wbasisCR)\
WfdParCRm<-fdPar(Wfd0CRm,1,1)\
registerlistCRm<-register.fd(accelmeanfdUNm,accelfdUNm,WfdParCRm)\
accelfdCRm<-registerlistCRm$regfd\
accelmeanfdCRm<-mean(accelfdCRm)

drv2sq_funm<-curvature_fda(heightbasis,accelmeanfdCRm)

OD_m<-ODsplines(n=n ,curvature=drv2sq_funm)

#### #OD for both genders using function 'ODsplines_sub'
OD_b<-ODsplines_sub(n=n, curvature1=drv2sq_funf,curvature2=drv2sq_funm)

#### #plot
plot(accelmeanfdCRf,xlim=c(0,1),xaxt="n",yaxt="n",
     main="",xlab="Age (years)",ylim=c(-1200,500))\
lines(accelmeanfdCRf,lty=1)\
lines(accelmeanfdCRm,lty=2)\
axis(1, at=c(seq(0,0.94,length.out =9),1) ,labels=c(1,3,5,7,9,11,13,15,17,18),las=1)\
axis(2, at=0 ,labels=0,las=1)

abline(h=0,col="gray")\
abline(h=-50,col="gray")\
abline(h=-100,col="gray")

points(OD_b,rep(-0,n),pch=8)\
points(OD_f,rep(-50,n),pch=15)\
points(OD_m,rep(-100,n),pch=16)\
legend(0.8,-900,legend=c("females","males"),lty=c(1,2),bty = "n")\
legend(0.85,-1000,legend=c("both","females","males"),pch=c(8,15:16),bty = "n")
