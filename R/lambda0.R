lambda0.matrix=function(x,lambda=0,maxit=100,trace.it=FALSE,thresh=1e-5){
  ismiss=is.na(x)
  if(any(ismiss))x[ismiss]=0
  svd(x)$d[1]
}

lambda0.Incomplete=function(x,lambda=0,maxit=100,trace.it=FALSE,thresh=1e-5){
  fit=svd.als(x,rank.max=2,lambda=lambda,maxit=maxit,trace.it=trace.it, thresh=thresh)
  lam0=fit$d[1]
  if(lam0<=thresh)warning("lambda too large; lambda0 is smaller than lambda")
  lam0+lambda
}

#' compute the smallest value for \code{lambda} such that
#' \code{softImpute(x,lambda)} returns the zero solution.
#' 
#' this determines the "starting" lambda for a sequence of values for
#' \code{softImpute}, and all nonzero solutions would require a smaller value
#' for \code{lambda}.
#' 
#' It is the largest singular value for the matrix, with zeros replacing
#' missing values. It uses \code{svd.als} with \code{rank=2}.
#' 
#' @aliases lambda0 lambda0,Incomplete-method lambda0,SparseplusLowRank-method
#' lambda0,sparseMatrix-method
#' @param x An m by n matrix. Large matrices can be in "sparseMatrix" format,
#' as well as "SparseplusLowRank". The latter arise after centering sparse
#' matrices, for example with \code{biScale}, as well as in applications such
#' as \code{softImpute}.
#' @param lambda As in \code{svd.als}, using a value for \code{lambda} can
#' speed up iterations. As long as the solution is not zero, the value returned
#' adds back this value.
#' @param maxit maximum number of iterations.
#' @param trace.it with \code{trace.it=TRUE}, convergence progress is reported.
#' @param thresh convergence threshold, measured as the relative changed in the
#' Frobenius norm between two successive estimates.
#' @return a single number, the largest singular value
#' @author Trevor Hastie, Rahul Mazumder\cr Maintainer: Trevor Hastie
#' \email{hastie@@stanford.edu}
#' @seealso \code{softImpute},\code{Incomplete}, and \code{svd.als}.
#' @references Rahul Mazumder, Trevor Hastie and Rob Tibshirani (2010)
#' \emph{Spectral Regularization Algorithms for Learning Large Incomplete
#' Matrices}, \url{https://hastie.su.domains/Papers/mazumder10a.pdf}\cr
#' \emph{ Journal of Machine Learning Research 11 (2010) 2287-2322}
#' @keywords models array multivariate
#' @examples
#' 
#' set.seed(101)
#' n=200
#' p=100
#' J=50
#' np=n*p
#' missfrac=0.3
#' x=matrix(rnorm(n*J),n,J)%*%matrix(rnorm(J*p),J,p)+matrix(rnorm(np),n,p)/5
#' ix=seq(np)
#' imiss=sample(ix,np*missfrac,replace=FALSE)
#' xna=x
#' xna[imiss]=NA
#' lambda0(xna)
#' @export
setGeneric("lambda0",lambda0.matrix)
setMethod("lambda0","Incomplete",lambda0.Incomplete)
setMethod("lambda0","SparseplusLowRank",lambda0.Incomplete)
setMethod("lambda0","sparseMatrix",lambda0.Incomplete)
