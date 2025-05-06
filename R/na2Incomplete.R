na2Incomplete=function(from){
    d=as.integer(dim(from))
  rows=rep(1:d[1],d[2])
  cols=rep(1:d[2],rep(d[1],d[2]))
  nisna=!is.na(from)  
  new("Incomplete",as(sparseMatrix(i=rows[nisna],j=cols[nisna],x=from[nisna],dims=d),"dgCMatrix"))
  }

#' @export
setAs("matrix","Incomplete",na2Incomplete)

sparse2Incomplete=function(from) new("Incomplete",as(from,"dgCMatrix"))

#' @export
setAs("sparseMatrix","Incomplete",sparse2Incomplete)

 
#' create a matrix of class \code{Incomplete}
#' 
#' creates an object of class \code{Incomplete}, which inherits from class
#' \code{dgCMatrix}, a specific instance of class \code{sparseMatrix}
#' 
#' The matrix is represented in sparse-matrix format, except the "zeros"
#' represent missing values. Real zeros are represented explicitly as values.
#' 
#' @aliases Incomplete coerce,matrix-method
#' @param i row indices
#' @param j column indices
#' @param x a vector of values
#' @return a matrix of class \code{Incomplete} which inherits from class
#' \code{dgCMatrix}
#' @author Trevor Hastie and Rahul Mazumder
#' @seealso \code{softImpute}
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
#' xnaC=as(xna,"Incomplete")
#' ### here we do it a different way to demonstrate Incomplete
#' ### In practise the observed values are stored in this market-matrix format.
#' i = row(xna)[-imiss]
#' j = col(xna)[-imiss]
#' xnaC=Incomplete(i,j,x=x[-imiss])
#' @export
Incomplete=function(i,j,x){
    new("Incomplete",as(sparseMatrix(i=i,j=j,x=x),"dgCMatrix"))
      }

  
