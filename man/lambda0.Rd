\name{lambda0}
\alias{lambda0}
\alias{lambda0,Incomplete-method}
\alias{lambda0,SparseplusLowRank-method}
\alias{lambda0,sparseMatrix-method}
\title{
compute the smallest value for \code{lambda} such that
\code{softImpute(x,lambda)} returns the zero solution.
}
\description{
this determines the "starting" lambda for a sequence of values for
\code{softImpute}, and all nonzero solutions would require a smaller
value for \code{lambda}.
}
\usage{
lambda0(x, lambda = 0, maxit = 100, trace.it = FALSE, thresh = 1e-05)
}
\arguments{
  \item{x}{
An m by n matrix. Large matrices can be in "sparseMatrix" format, as
well as "SparseplusLowRank". The latter arise after centering sparse
matrices, for example with \code{biScale}, as well as in applications
such as \code{softImpute}.
}
The remaining arguments only apply to matrices \code{x} in
\code{"sparseMatrix"}, \code{"Incomplete"}, or \code{"SparseplusLowRank"} format. 
\item{lambda}{
  As in \code{svd.als}, using a value for \code{lambda} can speed up
  iterations. As long as the solution is not zero, the value returned
  adds back this value.
  }
  \item{maxit}{
    maximum number of iterations.
  }
  \item{trace.it}{
with \code{trace.it=TRUE}, convergence progress is reported.
}
  \item{thresh}{
convergence threshold, measured as the relative changed in the Frobenius
norm between two successive estimates.
}
}
\details{
   It is the largest singular value for the matrix,
   with zeros replacing missing values. It uses \code{svd.als} with
   \code{rank=2}.
}
\value{
a single number, the largest singular value
}
\references{Rahul Mazumder, Trevor Hastie and Rob Tibshirani (2010)
  \emph{Spectral Regularization Algorithms for Learning Large Incomplete
Matrices},
  \url{https://web.stanford.edu/~hastie/Papers/mazumder10a.pdf}\cr
 \emph{ Journal of Machine Learning Research 11 (2010) 2287-2322}
}
\author{Trevor Hastie, Rahul Mazumder\cr
  Maintainer: Trevor Hastie  \email{hastie@stanford.edu}
}
\seealso{\code{softImpute},\code{Incomplete}, and \code{svd.als}.
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
lambda0(xna)
}
\keyword{models}
\keyword{array}
\keyword{multivariate}
