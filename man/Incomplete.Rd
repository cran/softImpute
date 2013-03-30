\name{Incomplete}
\alias{Incomplete}
\alias{coerce,matrix-method}
\title{
create a matrix of class \code{Incomplete}
}
\description{
creates an object of class \code{Incomplete}, which inherits from
class \code{dgCMatrix}, a specific instance of class \code{sparseMatrix}
}
\usage{
Incomplete(i, j, x)
}
\arguments{
  \item{i}{row indices
}
  \item{j}{column indices
}
  \item{x}{
a vector of values
}
}
\details{
The matrix is represented in sparse-matrix format, except the "zeros"
represent missing values. Real zeros are represented explicitly as values.
}
\value{
a matrix of class \code{Incomplete} which inherits from
class  \code{dgCMatrix} }
\author{
Trevor Hastie and Rahul Mazumder
}
\seealso{\code{softImpute}
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
xnaC=as(xna,"Incomplete")
### here we do it a different way to demonstrate Incomplete
### In practise the observed values are stored in this market-matrix format.
i = row(xna)[-imiss]
j = col(xna)[-imiss]
xnaC=Incomplete(i,j,x=x[-imiss])
  }
\keyword{models}
\keyword{array}
\keyword{multivariate}
