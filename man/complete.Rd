\name{complete}
\alias{complete}
\alias{impute}
\alias{complete,Incomplete-method}
\alias{complete,matrix-method}
\title{
make predictions from an svd object
}
\description{
These functions produce predictions from the low-rank solution of \code{softImpute}
}
\usage{
complete(x, object, unscale = TRUE)
impute(object, i, j, unscale = TRUE)
}
\arguments{
  \item{x}{
a matrix with NAs or a matrix of class \code{"Incomplete"}.
}
  \item{object}{
an svd object with components u, d and v
}
 \item{i}{vector of row indices for the locations to be predicted}
 \item{j}{vector of column indices for the locations to be predicted}
  \item{unscale}{
if \code{object} has \code{biScale} attributes, and \code{unscale=TRUE},
the imputations reversed the centering and scaling on the predictions.
}
}
\details{
\code{impute} returns a vector of predictions, using the reconstructed
low-rank matrix representation represented by \code{object}. It is used by complete,
which returns a complete matrix with all the missing values imputed.
}
\value{
Either a vector of predictions or a complete matrix. WARNING: if
\code{x} has large dimensions, the matrix returned by \code{complete}
might be too large.}
\author{
Trevor Hastie
}
\seealso{\code{softImpute}, \code{biScale} and \code{Incomplete}
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
complete(xna,fit1)
}
\keyword{models}
\keyword{array}
\keyword{multivariate}
