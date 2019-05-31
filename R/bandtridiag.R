#' Generate band matrix with three diagonals on the lower side
#'
#'Generate band matrix of dimention p*(p-2) that has nonzero elements only on the main diagonal, the first two diagonals below the main diagonal.
#'
#' @param main  vector of length p-1, the elements in the main diagonal.
#' @param middle vector of length p-1, the elements in the first diagonal below the main diagonal.
#' @param lower vector of length p-1, the elements in the second diagonal below the main diagonal.
#' @examples lower<-1:3
#' middle<-4:6
#' main<-7:9
#' bandtridiag(lower, middle, main)
#'
#' @export
bandtridiag<-function(lower, middle, main){
  out <- matrix(0,length(main),length(main)+2)
  diag(out) <- main
  indx <- seq.int(length(lower))
  out[cbind(indx,indx+1)] <- middle
  out[cbind(indx,indx+2)] <- lower
  return(t(out))
}
