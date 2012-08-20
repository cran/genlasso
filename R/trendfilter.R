# We compute the solution path of the trend filtering problem:
#
# \hat{\beta}(\lambda) =
# \argmin_\beta \|y - X \beta|_2^2 + \lambda\|D \beta\|_1,
#
# where D is (p-k-1) x p is the discrete difference operator of
# order k+1, and X is n x p and full column rank. The solution is
# a piecewise polynomial of degree k, with adaptively chosen knots.

trendfilter <- function(y, X, z, ord=1, approx=FALSE, maxsteps=2000,
                        minlam=0, tol=1e-11, verbose=FALSE) {
  cl = match.call()

  if (missing(y)) stop("y is missing.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (missing(z)) z = NULL
  if (!is.null(z) && !is.numeric(z)) stop("z must be numeric.")
  if (!is.null(z) && is.unsorted(z)) stop("z must be in increasing order.")
  if (missing(X)) X = NULL
  if (!is.null(X) && !is.matrix(X)) stop("X must be a matrix.")
  if (!is.null(z) && !is.null(X)) stop("z cannot be used with a design matrix X.")
  if (ord<0 || round(ord)!=ord) stop("ord must be a nonnegative integer.")
  if (length(y) <= ord+1) stop("Not enough data points to fit a trend of the given order [need length(y) > ord+1].")
  if (ord>3) warning(paste("For numerical stability, it is not recommended to run",
                           "trend filtering with a polynomial order larger than 3."))
  
  if (is.null(X)) {
    n = length(y)    
    if (is.null(z)) D = getDtfSparse(n,ord)
    else D = getDtfPosSparse(n,ord,z)
    out = dualpathWideSparse(y,D,NULL,approx,maxsteps,minlam,tol,verbose)
  }

  else {
    if (length(y)!=nrow(X)) stop("Dimensions don't match [length(y) != nrow(X)].")
    p = ncol(X)
    D = getDtfSparse(p,ord)
    out = dualpathTrendX(y,X,D,ord,approx,maxsteps,minlam,tol,verbose)
  }
  
  out$trendorder = ord
  if (!is.null(z)) out$z = z
  out$call = cl
  
  class(out) = c("trendfilter", "genlasso", "list")
  return(out)
}
