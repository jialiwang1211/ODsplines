% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ODsplines_sub.R
\name{ODsplines_sub}
\alias{ODsplines_sub}
\title{Optimal design points for two populations}
\usage{
ODsplines_sub(n, eta = 1, curvature1, curvature2, interv = rep(0.001, n
  - 1), seed = 1)
}
\arguments{
\item{n}{number of design points.}

\item{eta}{smoothing parameter.}

\item{curvature1}{square of prior curvature of population 1, function of t.}

\item{curvature2}{square of prior curvature of population 2, function of t.}

\item{interv}{a vector of length n-1, the mininal intervals between consecutive design points.}

\item{seed}{an integer.}
}
\description{
Derive the optimal design points of the adaptive smoothing splines with different prior curvatures from two subpopulation.
}
\examples{
n<-7
curvature1<-curvature_logistic(theta1=1,theta2=-12,theta3=5)
curvature2<-curvature_logistic(theta1=1,theta2=-9,theta3=6)
OD_t<-ODsplines_sub(n=n,curvature1=curvature1,curvature2=curvature2)

}
