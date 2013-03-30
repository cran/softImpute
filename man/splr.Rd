\name{splr}
\alias{splr}
\title{
create a \code{SparseplusLowRank} object
}
\description{
create an object of class  \code{SparseplusLowRank} which can be
efficiently stored and for which efficient linear algebra operations are possible.
}
\usage{
splr(x, a = NULL, b = NULL)
}
\arguments{
  \item{x}{
sparse matrix with dimension say m x n
}
  \item{a}{
matrix with m rows and number of columns  r less than \code{min(dim(x))}
}
  \item{b}{
matrix with n rows and number of columns  r less than \code{min(dim(x))}
}
}
\value{
an object of S4 class \code{SparseplusLowRank} is returned with slots
\code{x}, \code{a} and \code{b}
}
\author{
Trevor Hastie
}
\seealso{\code{SparseplusLowRank-class}, \code{softImpute}
}
\examples{
x=matrix(sample(c(3,0),15,replace=TRUE),5,3)
x=as(x,"sparseMatrix")
a=matrix(rnorm(10),5,2)
b=matrix(rnorm(6),3,2)
new("SparseplusLowRank",x=x,a=a,b=b)
splr(x,a,b)
}
\keyword{classes}
\keyword{models}
\keyword{multivariate}