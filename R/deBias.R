#' Recompute the \code{$d} component of a \code{"softImpute"} object through
#' regression.
#' 
#' \code{softImpute} uses shrinkage when completing a matrix with missing
#' values. This function debiases the singular values using ordinary least
#' squares.
#' 
#' Treating the \code{"d"} values as parameters, this function recomputes them
#' by linear regression.
#' 
#' @param x matrix with missing entries, or a matrix of class
#' \code{"Incomplete"}
#' @param svdObject an SVD object, the output of \code{softImpute}
#' @return An svd object is returned, with components "u", "d", and "v".
#' @author Trevor Hastie\cr Maintainer: Trevor Hastie
#' \email{hastie@@stanford.edu}
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
#' fit1=softImpute(xna,rank=50,lambda=30)
#' fit1d=deBias(xna,fit1)
#' @export
deBias=function(x,svdObject){
  if(!inherits(x,"Incomplete"))x=as(x,"Incomplete")
  irow=x@i
  pcol=x@p
  x=x@x
  u=as.matrix(svdObject$u)
  v=as.matrix(svdObject$v)

  dd=dim(u)
  nnrow=as.integer(dd[1])
  nncol=as.integer(nrow(v))
  nrank=dd[2]
  storage.mode(u)="double"
  storage.mode(v)="double"
  storage.mode(irow)="integer"
  storage.mode(pcol)="integer"
  storage.mode(x)="double"
  nomega=as.integer(length(irow))
  out=.Fortran("plusregC",
           nnrow,nncol,nrank,u,v,irow,pcol,nomega,x,
           zy=double(nrank),zz=matrix(double(nrank*nrank),nrank,nrank),
           PACKAGE="softImpute"
           )
  zz=out$zz
  which=lower.tri(zz)
  zz[which]=t(zz)[which]
  zy=out$zy
  d=solve(zz,zy)
  sd=sign(d)
  if(any(sd<0)){
    d=abs(d)
    o=rev(order(d))
    v=scale(v,FALSE,sd)[,o]
    u=u[,o]
    d=d[o]
  }
  list(u=u,v=v,d=d)
}

