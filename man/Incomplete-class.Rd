\name{Incomplete-class}
\Rdversion{1.1}
\docType{class}
\alias{Incomplete-class}
\alias{as.matrix,Incomplete-method}
\alias{coerce,matrix,Incomplete-method}
\alias{coerce,sparseMatrix,Incomplete-method}

\title{Class \code{"Incomplete"}}
\description{
a sparse matrix inheriting from class \code{dgCMatrix} with the NAs
represented as zeros
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Incomplete", ...)}.
or by calling the function \code{Incomplete}
}
\section{Slots}{
  \describe{
    \item{\code{i}:}{Object of class \code{"integer"} ~~ }
    \item{\code{p}:}{Object of class \code{"integer"} ~~ }
    \item{\code{Dim}:}{Object of class \code{"integer"} ~~ }
    \item{\code{Dimnames}:}{Object of class \code{"list"} ~~ }
    \item{\code{x}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{factors}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{dgCMatrix}"}, directly.
Class \code{"\linkS4class{CsparseMatrix}"}, by class "dgCMatrix", distance 2.
Class \code{"\linkS4class{dsparseMatrix}"}, by class "dgCMatrix", distance 2.
Class \code{"\linkS4class{generalMatrix}"}, by class "dgCMatrix", distance 2.
Class \code{"\linkS4class{dMatrix}"}, by class "dgCMatrix", distance 3.
Class \code{"\linkS4class{sparseMatrix}"}, by class "dgCMatrix", distance 3.
Class \code{"\linkS4class{compMatrix}"}, by class "dgCMatrix", distance 3.
Class \code{"\linkS4class{Matrix}"}, by class "dgCMatrix", distance 4.
}
\section{Methods}{
  \describe{
    \item{as.matrix}{\code{signature(x = "Incomplete")}: ... }
    \item{coerce}{\code{signature(from = "matrix", to = "Incomplete")}: ... }
    \item{complete}{\code{signature(x = "Incomplete")}: ... }
	 }
}
\author{
Trevor Hastie and Rahul Mazumder
}
\seealso{\code{biScale},\code{softImpute},\code{Incomplete},\code{impute},\code{complete}
}
\examples{
showClass("Incomplete")
set.seed(101)
n=200
p=100
J=50
np=n*p
missfrac=0.3
x=matrix(rnorm(n*J),n,J)\%*\%matrix(rnorm(J*p),J,p)+matrix(rnorm(np),n,p)/5
ix=seq(np)
imiss=sample(ix,np*missfrac,replace=FALSE)
xna=x
xna[imiss]=NA
xnaC=as(xna,"Incomplete")
}
\keyword{classes}
