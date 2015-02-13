\name{SparseplusLowRank-class}
\Rdversion{1.1}
\docType{class}
\alias{SparseplusLowRank-class}
\alias{\%*\%,ANY,SparseplusLowRank-method}
\alias{\%*\%,SparseplusLowRank,ANY-method}
\alias{\%*\%,Matrix,SparseplusLowRank-method}
\alias{\%*\%,SparseplusLowRank,Matrix-method}
\alias{as.matrix,SparseplusLowRank-method}
\alias{colMeans,SparseplusLowRank-method}
\alias{colSums,SparseplusLowRank-method}
\alias{dim,SparseplusLowRank-method}
\alias{norm,SparseplusLowRank,character-method}
\alias{rowMeans,SparseplusLowRank-method}
\alias{rowSums,SparseplusLowRank-method}

\title{Class \code{"SparseplusLowRank"}}
\description{
A structured matrix made up of a sparse part plus a low-rank part, all
which can be stored and operated on efficiently.
}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{new("SparseplusLowRank", ...)}
or by a call to \code{splr}
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"sparseMatrix"}}
    \item{\code{a}:}{Object of class \code{"matrix"}}
    \item{\code{b}:}{Object of class \code{"matrix"}}
  }
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "ANY", y = "SparseplusLowRank")}: ... }
    \item{\%*\%}{\code{signature(x = "SparseplusLowRank", y = "ANY")}: ... }
    \item{\%*\%}{\code{signature(x = "Matrix", y = "SparseplusLowRank")}: ... }
    \item{\%*\%}{\code{signature(x = "SparseplusLowRank", y = "Matrix")}: ... }
    \item{as.matrix}{\code{signature(x = "SparseplusLowRank")}: ... }
    \item{colMeans}{\code{signature(x = "SparseplusLowRank")}: ... }
    \item{colSums}{\code{signature(x = "SparseplusLowRank")}: ... }
    \item{dim}{\code{signature(x = "SparseplusLowRank")}: ... }
    \item{norm}{\code{signature(x = "SparseplusLowRank", type = "character")}: ... }
    \item{rowMeans}{\code{signature(x = "SparseplusLowRank")}: ... }
    \item{rowSums}{\code{signature(x = "SparseplusLowRank")}: ... }
    \item{svd.als}{\code{signature(x = "SparseplusLowRank")}: ... }
	 }
}
\author{
Trevor Hastie and Rahul Mazumder
}
\seealso{\code{softImpute},\code{splr}
}
\examples{
showClass("SparseplusLowRank")
x=matrix(sample(c(3,0),15,replace=TRUE),5,3)
x=as(x,"sparseMatrix")
a=matrix(rnorm(10),5,2)
b=matrix(rnorm(6),3,2)
new("SparseplusLowRank",x=x,a=a,b=b)
splr(x,a,b)
}
\keyword{classes}
