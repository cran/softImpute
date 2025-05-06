#' Softimpute for matrix completion
#'
#' SoftImpute solves the following problem for a matrix \eqn{X} with missing
#' entries: \eqn{\min||X-M||_o^2 +\lambda||M||_*.} Here \eqn{||\cdot||_o} is
#' the Frobenius norm, restricted to the entries corresponding to the
#' non-missing entries of \eqn{X}, and \eqn{||M||_*} is the nuclear norm of
#' \eqn{M} (sum of singular values).  For full details of the "svd" algorithm
#' are described in the reference below.  The "als" algorithm will be described
#' in a forthcoming article. Both methods employ special sparse-matrix tricks
#' for large matrices with many missing values. This package creates a new
#' sparse-matrix class \code{"SparseplusLowRank"} for matrices of the form
#' \deqn{x+ab',} where \eqn{x} is sparse and \eqn{a} and \eqn{b} are tall
#' skinny matrices, hence \eqn{ab'} is low rank. Methods for efficient left and
#' right matrix multiplication are provided for this class. For large matrices,
#' the function \code{Incomplete()} can be used to build the appropriate sparse
#' input matrix from market-format data.
#' @references Rahul Mazumder, Trevor Hastie and Rob Tibshirani (2010)
#' \emph{Spectral Regularization Algorithms for Learning Large Incomplete
#' Matrices}, \url{https://hastie.su.domains/Papers/mazumder10a.pdf}
#' \emph{ Journal of Machine Learning Research 11 (2010) 2287-2322}
#'
#' @name softImpute-package
#' @import methods
#' @import Matrix
#' @importFrom utils packageDescription
#' @importFrom stats rnorm weighted.mean
#' @useDynLib softImpute
#' @author Trevor Hastie and Rahul Mazumder
#'
#' Maintainer: Trevor Hastie<hastie@stanford.edu>
#' @keywords internal
"_PACKAGE"



