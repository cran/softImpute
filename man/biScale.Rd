\name{biScale}
\alias{biScale}
\title{
standardize a matrix to have optionally row means zero and variances one,
and/or column means zero and variances one.
}
\description{
A function for standardizing a matrix in a symmetric
fashion. Generalizes the \code{scale} function in R. Works with
matrices with NAs, matrices of class "Incomplete", and matrix in
"sparseMatrix" format.
}
\usage{
biScale(x, maxit = 20, thresh = 1e-09, row.center = TRUE, row.scale =TRUE,
        col.center = TRUE, col.scale = TRUE, trace = FALSE)
}
\arguments{
  \item{x}{
matrix, possibly with NAs, also of class "Incomplete" or "sparseMatrix" format.
}
  \item{maxit}{
When both row and column centering/scaling is requested, iteration is
may be necessary
}
  \item{thresh}{
Convergence threshold
}
  \item{row.center}{
if \code{row.center==TRUE} (the default), row centering will be performed resulting in a matrix
with row means zero. If \code{row.center} is a vector, it will be used
to center the rows. If \code{row.center=FALSE} nothing is done. See details for
more info.
}
  \item{row.scale}{
if \code{row.scale==TRUE}, the rows are scaled (after possibly
centering, to have variance one. Alternatively, if a positive vector is
supplied, it is used for row centering. See details for
more info.
}
  \item{col.center}{
Similar to \code{row.center}
}
  \item{col.scale}{
Similar to \code{row.scale}
}
  \item{trace}{
with \code{trace=TRUE}, convergence progress is reported, when iteration
is needed.
}
}
\details{
  This function computes a transformation
  \deqn{\frac{X_{ij}-\alpha_i-\beta_j}{\gamma_i\tau_j}}
  to transform the matrix \eqn{X}. It uses an iterative algorithm based
  on "method-of-moments". At each step, all but one of the parameter
  vectors is fixed, and the remaining vector is computed to solve the
  required condition. Although in genereal this is not guaranteed to
  converge,
  it mostly does, and quite rapidly. When there are convergence
  problems, remove some of the required constraints. When any of the
  row/column centers or scales are provided, they are used rather than
  estimated in the above model.
}
\value{
  A matrix like \code{x} is returned, with attributes:
  \item{biScale:row}{a list with elements \code{"center"} and
    \code{"scale"} (the \eqn{alpha} and \eqn{gamma} above. If no
    centering was done, the center component will be a vector of
    zeros. Likewise, of no row scaling was done, the scale component
    will be a vector of ones.}
  \item{biScale:column}{Same details as \code{biScale:row}}
  For matrices with missing values, the constraints apply to the
  non-missing entries. If \code{x} is of class \code{"sparseMatrix"},
  then the sparsity is maintained, and an object of class
  \code{"SparseplusLowRank"} is returned, such that the low-rank part
  does the centering.
  }
\author{Trevor Hastie, with help from Andreas Buja and Steven Boyd\cr,
  Maintainer: Trevor Hastie  \email{hastie@stanford.edu}
}
\note{This function will be described in detail in a forthcoming paper}
\seealso{\code{softImpute},\code{Incomplete},\code{lambda0},\code{impute},\code{complete},
  and class \code{"SparseplusLowRank"}
    }

\examples{
set.seed(101)
n=200
p=100
J=50
np=n*p
missfrac=0.3
x=matrix(rnorm(n*J),n,J)\%*\%matrix(rnorm(J*p),J,p)+matrix(rnorm(np),n,p)/5
xc=biScale(x)
ix=seq(np)
imiss=sample(ix,np*missfrac,replace=FALSE)
xna=x
xna[imiss]=NA
xnab=biScale(xna,row.scale=FALSE,trace=TRUE)
xnaC=as(xna,"Incomplete")
xnaCb=biScale(xnaC)
nnz=trunc(np*.3)
inz=sample(seq(np),nnz,replace=FALSE)
i=row(x)[inz]
j=col(x)[inz]
x=rnorm(nnz)
xS=sparseMatrix(x=x,i=i,j=j)
xSb=biScale(xS)
class(xSb)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{array}
\keyword{multivariate}
