#' compute a low rank soft-thresholded svd by alternating orthogonal ridge
#' regression
#' 
#' fit a low-rank svd to a complete matrix by alternating orthogonal ridge
#' regression. Special sparse-matrix classes available for very large matrices,
#' including "SparseplusLowRank" versions for row and column centered sparse
#' matrices.
#' 
#' This algorithm solves the problem \deqn{\min ||X-M||_F^2 +\lambda ||M||_*}
#' subject to \eqn{rank(M)\leq r}, where \eqn{||M||_*} is the nuclear norm of
#' \eqn{M} (sum of singular values). It achieves this by solving the related
#' problem \deqn{\min ||X-AB'||_F^2 +\lambda/2 (||A||_F^2+||B||_F^2)} subject
#' to \eqn{rank(A)=rank(B)\leq r}. The solution is a rank-restricted,
#' soft-thresholded SVD of \eqn{X}.
#' 
#' @aliases svd.als svd.als,sparseMatrix-method
#' svd.als,SparseplusLowRank-method
#' @param x An m by n matrix. Large matrices can be in "sparseMatrix" format,
#' as well as "SparseplusLowRank". The latter arise after centering sparse
#' matrices, for example with \code{biScale}, as well as in applications such
#' as \code{softImpute}.
#' @param rank.max The maximum rank for the solution. This is also the
#' dimension of the left and right matrices of orthogonal singular vectors.
#' 'rank.max' should be no bigger than 'min(dim(x)'.
#' @param lambda The regularization parameter. \code{lambda=0} corresponds to
#' an accelerated version of the orthogonal QR-algorithm. With \code{lambda>0}
#' the algorithm amounts to alternating orthogonal ridge regression.
#' @param thresh convergence threshold, measured as the relative changed in the
#' Frobenius norm between two successive estimates.
#' @param maxit maximum number of iterations.
#' @param trace.it with \code{trace.it=TRUE}, convergence progress is reported.
#' @param warm.start an svd object can be supplied as a warm start. If the
#' solution requested has higher rank than the warm start, the additional
#' subspace is initialized with random Gaussians (and then orthogonalized wrt
#' the rest).
#' @param final.svd Although in theory, this algorithm converges to the
#' solution to a nuclear-norm regularized low-rank matrix approximation
#' problem, with potentially some singular values equal to zero, in practice
#' only near-zeros are achieved. This final step does one more iteration with
#' \code{lambda=0}, followed by soft-thresholding.
#' @return An svd object is returned, with components "u", "d", and "v".
#' \item{u}{an m by \code{rank.max} matrix with the left orthogonal singular
#' vectors} \item{d}{a vector of length \code{rank.max} of soft-thresholded
#' singular values} \item{v}{an n by \code{rank.max} matrix with the right
#' orthogonal singular vectors}
#' @author Trevor Hastie, Rahul Mazumder\cr Maintainer: Trevor Hastie
#' \email{hastie@@stanford.edu}
#' @seealso \code{biScale}, \code{softImpute}, \code{Incomplete},
#' \code{lambda0}, \code{impute}, \code{complete}
#' @references Rahul Mazumder, Trevor Hastie and Rob Tibshirani (2010)
#' \emph{Spectral Regularization Algorithms for Learning Large Incomplete
#' Matrices}, \url{https://hastie.su.domains/Papers/mazumder10a.pdf}\cr
#' \emph{ Journal of Machine Learning Research 11 (2010) 2287-2322}
#' @keywords models array multivariate
#' @examples
#' 
#' #create a matrix and run the algorithm
#' set.seed(101)
#' n=100
#' p=50
#' J=25
#' np=n*p
#' x=matrix(rnorm(n*J),n,J)%*%matrix(rnorm(J*p),J,p)+matrix(rnorm(np),n,p)/5
#' fit=svd.als(x,rank=25,lambda=50)
#' fit$d
#' pmax(svd(x)$d-50,0)
#' # now create a sparse matrix and do the same
#' nnz=trunc(np*.3)
#' inz=sample(seq(np),nnz,replace=FALSE)
#' i=row(x)[inz]
#' j=col(x)[inz]
#' x=rnorm(nnz)
#' xS=sparseMatrix(x=x,i=i,j=j)
#' fit2=svd.als(xS,rank=20,lambda=7)
#' fit2$d
#' pmax(svd(as.matrix(xS))$d-7,0)
#' @export
svd.als=function(x, rank.max=2, lambda=0,thresh = 1e-05, maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE){
  if(rank.max>(rmax<-min(dim(x)))){
    rank.max=rmax
    warning(paste("rank.max should not exceed min(dim(x)); changed to ",rmax))
  }
  ismiss=is.na(x)
  if(any(ismiss))stop("NAs in x; use softImpute instead")
  this.call=match.call()
  out=simpute.als(x,J=rank.max,thresh,lambda,maxit,trace.it,warm.start,final.svd)
  attr(out,"call")=this.call
  attr(out,"lambda")=lambda
  out
}

svd.als.sparse=function(x, rank.max=2, lambda=0,thresh = 1e-05, maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE){
   this.call=match.call()
  out=Ssvd.als(x,J=rank.max,thresh,lambda,maxit,trace.it,warm.start,final.svd)
  attr(out,"call")=this.call
  attr(out,"lambda")=lambda
  out
}

  setGeneric("svd.als",svd.als)
  setMethod("svd.als","sparseMatrix",svd.als.sparse)
  setMethod("svd.als","SparseplusLowRank",svd.als.sparse)
  
