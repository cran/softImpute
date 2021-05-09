\name{softImpute}
\alias{softImpute}
\title{impute missing values for a matrix via nuclear-norm regularization.
}
\description{fit a low-rank matrix approximation to a matrix with
  missing values via nuclear-norm regularization. The algorithm works
  like EM, filling in the missing values with the current guess, and
  then solving the optimization problem on the complete matrix using a
  soft-thresholded SVD. Special sparse-matrix classes available for very
  large matrices.
}
\usage{
softImpute(x, rank.max = 2, lambda = 0, type = c("als", "svd"), thresh = 1e-05,
           maxit = 100, trace.it = FALSE, warm.start = NULL, final.svd = TRUE)
}
\arguments{
  \item{x}{An m by n matrix with NAs. For large matrices can be of class
    \code{"Incomplete"}, in which case the missing values are
    represented as pseudo zeros leading to dramatic storage
    reduction. \code{x} can have been centered and scaled via
    \code{biScale}, and this information is carried along with the solution.
}
  \item{rank.max}{
This restricts the rank of the solution. If sufficiently large, and with
\code{type="svd"}, the solution solves the nuclear-norm convex
matrix-completion problem. In this case the number of nonzero singular
values returned will be less than or equal to \code{rank.max}. If smaller
ranks are used, the solution is not guaranteed to solve the problem,
although still results in good local minima. \code{rank.max} should be
no bigger than \code{min(dim(x)-1}.
}
  \item{lambda}{
nuclear-norm regularization parameter. If \code{lambda=0}, the algorithm
reverts to "hardImpute", for which convergence is typically slower, and
to local minimum. Ideally \code{lambda} should be chosen so that the solution
reached has rank slightly less than \code{rank.max}. See also
\code{lambda0()} for computing the smallest \code{lambda} with a zero solution.
}
  \item{type}{two algorithms are implements, \code{type="svd"} or
  the default  \code{type="als"}. The "svd" algorithm repeatedly computes
  the svd of the completed matrix, and soft thresholds its singular
  values. Each new soft-thresholded svd is used to re-impute the missing
  entries. For large matrices of class \code{"Incomplete"}, the svd is
  achieved by an efficient form of alternating orthogonal ridge
  regression. The "als" algorithm uses this same alternating ridge
  regression, but updates the imputation at each step, leading to quite
  substantial speedups in some cases. The "als" approach does not
  currently have the same theoretical  convergence guarantees as the
  "svd" approach.
  }
  \item{thresh}{
convergence threshold, measured as the relative change in the Frobenius
norm between two successive estimates.
}
  \item{maxit}{
    maximum number of iterations.
  }
  \item{trace.it}{
with \code{trace.it=TRUE}, convergence progress is reported.
}
  \item{warm.start}{
an svd object can be supplied as a warm start. This is particularly
useful when constructing a path of solutions with decreasing values of
\code{lambda} and increasing \code{rank.max}. The previous solution can
be provided directly as a warm start for the next.
}
  \item{final.svd}{
only applicable to \code{type="als"}. The alternating ridge-regressions
do not lead to exact zeros. With the default \code{final.svd=TRUE}, at
the final iteration, a one step unregularized iteration is performed,
followed by soft-thresholding of the singular values, leading to hard zeros.
}
}
\details{
  SoftImpute solves the following problem for a matrix \eqn{X} with
  missing entries:
  \deqn{\min||X-M||_o^2 +\lambda||M||_*.}
Here \eqn{||\cdot||_o} is the Frobenius norm, restricted to the entries
corresponding to the
non-missing entries of \eqn{X}, and  \eqn{||M||_*} is the nuclear norm
of \eqn{M} (sum of singular values).  
  For full details of the "svd" algorithm are described in the reference
below.  The "als" algorithm will be described in a forthcoming
article. Both methods employ special sparse-matrix tricks for large
matrices with many missing values. This package creates a new
sparse-matrix class \code{"SparseplusLowRank"} for matrices of the form
\deqn{x+ab',} where \eqn{x} is sparse and \eqn{a} and \eqn{b} are tall
skinny matrices, hence \eqn{ab'} is low rank. Methods for efficient left
and right matrix multiplication are provided for this class. For large
matrices, the function \code{Incomplete()} can be used to build the
appropriate
sparse input matrix from market-format data.
}
\value{ An svd object is returned, with components "u", "d", and "v".
  If the solution has zeros in "d", the solution is truncated to rank one
  more than the number of zeros (so the zero is visible). If the input
  matrix had been centered and scaled by \code{biScale}, the scaling
  details are assigned as attributes inherited from the input matrix.
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
\seealso{\code{biScale}, \code{svd.als},\code{Incomplete},
  \code{lambda0}, \code{impute}, \code{complete}
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
###uses regular matrix method for matrices with NAs
fit1=softImpute(xna,rank=50,lambda=30)
###uses sparse matrix method for matrices of class "Incomplete"
xnaC=as(xna,"Incomplete")
fit2=softImpute(xnaC,rank=50,lambda=30)
###uses "svd" algorithm
fit3=softImpute(xnaC,rank=50,lambda=30,type="svd")
ximp=complete(xna,fit1)
### first scale xna
xnas=biScale(xna)
fit4=softImpute(xnas,rank=50,lambda=10)
ximp=complete(xna,fit4)
impute(fit4,i=c(1,3,7),j=c(2,5,10))
impute(fit4,i=c(1,3,7),j=c(2,5,10),unscale=FALSE)#ignore scaling and centering
  }
\keyword{models}
\keyword{array}
\keyword{multivariate}

