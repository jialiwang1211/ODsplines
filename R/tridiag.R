#' Generate tridiagonal matrix
#'
#'Generate square tridiagonal matrix of dimention p that has nonzero elements only on the main diagonal, the first diagonal below and above the main diagonal.
#'
#' @param main  vector of length p,the elements in the main diagonal.
#' @param upper vector of length p-1, the elements in the first diagonal above the main diagonal.
#' @param lower vector of length p-1, the elements in the first diagonal below the main diagonal.
#'
#' @examples main<-1:4
#' upper<-5:7
#' lower<-8:10
#' tridiag(upper, lower, main)
#'
#' @export
tridiag <- function(upper, lower, main){
  out <- matrix(0,length(main),length(main))
  diag(out) <- main
  indx <- seq.int(length(upper))
  out[cbind(indx+1,indx)] <- lower
  out[cbind(indx,indx+1)] <- upper
  return(out)
}
