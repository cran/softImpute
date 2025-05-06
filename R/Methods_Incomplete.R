#' Class \code{"Incomplete"}
#' 
#' a sparse matrix inheriting from class \code{dgCMatrix} with the NAs
#' represented as zeros
#' 
#' 
#' @name Incomplete-class
#' @aliases Incomplete-class as.matrix,Incomplete-method
#' coerce,matrix,Incomplete-method coerce,sparseMatrix,Incomplete-method
#' @docType class
#'
#' @section Objects from the Class: Objects can be created by calls of the form
#' `new("Incomplete", ...)` or by calling the function `Incomplete`
#'
#' @section Slots:
#' \describe{
#' \item{i}{Object of class `"integer"`}
#' \item{p}{Object of class `"integer"`}
#' \item{Dim}{Object of class `"integer"`}
#' \item{Dimnames}{Object of class `"list"`}
#' \item{x}{Object of class `"numeric"`}
#' \item{factors}{Object of class `"list"`}
#' }
#'
#' @section Methods:
#' \describe{
#' \item{as.matrix}{\code{signature(x = "Incomplete")}: ... }
#' \item{coerce}{\code{signature(from = "matrix", to = "Incomplete")}: ... }
#' \item{complete}{\code{signature(x = "Incomplete")}:... }
#' }
#' @author Trevor Hastie and Rahul Mazumder
#' @seealso
#' \code{biScale},\code{softImpute},\code{Incomplete},\code{impute},\code{complete}
#' @keywords classes
#' @examples
#' 
#' showClass("Incomplete")
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
#' @export
setClass("Incomplete","dgCMatrix")
as.matrix.Incomplete=function(x,...){
  x=as(x,"dgTMatrix")
  i=x@i
  j=x@j
  out=as.matrix(x)
  out[]=NA
  nrow=dim(x)[1]
  out[i+1 +nrow*j]=x@x
  out
}
setMethod("as.matrix","Incomplete",as.matrix.Incomplete)
complete.matrix=function(x,object,unscale=TRUE){
  nas=is.na(x)
  if(!any(nas))return(x)
  i=row(x)[nas]
  j=col(x)[nas]
  x[nas]=impute(object,i,j,unscale=unscale)
  x
}
scomplete = function(x,object,unscale=TRUE){
###x is class Incomplete
    x=as.matrix(x)##Fills in NAs where the zeros are
    complete(x,object,unscale=unscale)
  }

#' @export
setGeneric("complete",complete.matrix)

#' @export
setMethod("complete","Incomplete",scomplete)
