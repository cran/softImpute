#' make predictions from an svd object
#' 
#' These functions produce predictions from the low-rank solution of
#' \code{softImpute}
#' 
#' \code{impute} returns a vector of predictions, using the reconstructed
#' low-rank matrix representation represented by \code{object}. It is used by
#' complete, which returns a complete matrix with all the missing values
#' imputed.
#' 
#' @aliases complete impute complete,Incomplete-method complete,matrix-method
#' @param object an svd object with components u, d and v
#' @param i vector of row indices for the locations to be predicted
#' @param j vector of column indices for the locations to be predicted
#' @param unscale if \code{object} has \code{biScale} attributes, and
#' \code{unscale=TRUE}, the imputations reversed the centering and scaling on
#' the predictions.
#' @return Either a vector of predictions or a complete matrix. WARNING: if
#' \code{x} has large dimensions, the matrix returned by \code{complete} might
#' be too large.
#' @author Trevor Hastie
#' @seealso \code{softImpute}, \code{biScale} and \code{Incomplete}
#' @keywords models array multivariate
#' @examples
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
#' fit1=softImpute(xna,rank=50,lambda=30)
#' complete(xna,fit1)
#' @export
impute=function(object,i,j,unscale=TRUE){
  ###obj is a svd object produced by our Softimpute algorithms
  v=as.matrix(object$v)
 vd=v*outer(rep(1,nrow(v)),object$d)
 out= suv(as.matrix(object$u),vd,i,j)
 if(unscale){
   biats=attributes(object)
   if(any(grep("biScale",names(biats)))){
     out=out*biats[["biScale:row"]]$scale[i]*biats[["biScale:column"]]$scale[j]
     out=out+biats[["biScale:row"]]$center[i]+biats[["biScale:column"]]$center[j]
   }
 }
 out
}
