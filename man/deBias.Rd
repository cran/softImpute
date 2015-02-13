\name{deBias}
\alias{deBias}
\title{Recompute the \code{$d} component of a \code{"softImpute"} object
  through regression.
}
\description{\code{softImpute} uses shrinkage when completing a matrix with
  missing values. This function debiases the singular values using
  ordinary least squares.
}
\usage{
deBias(x, svdObject)
}
\arguments{
  \item{x}{matrix with missing entries, or a matrix of class \code{"Incomplete"}
}
  \item{svdObject}{an SVD object, the output of \code{softImpute}
}
}
\details{
Treating the \code{"d"} values as parameters, this function recomputes
them by linear regression.
}
\value{
   An svd object is returned, with components "u", "d", and "v".
  }
\author{Trevor Hastie\cr
  Maintainer: Trevor Hastie  \email{hastie@stanford.edu}
}
\examples{
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
fit1=softImpute(xna,rank=50,lambda=30)
fit1d=deBias(xna,fit1)
}
\keyword{models}
\keyword{array}
\keyword{multivariate}

