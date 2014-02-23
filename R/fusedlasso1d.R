# We compute the solution path of the trend filtering problem:
#
# \hat{\beta}(\lambda) =
# \argmin_\beta \|y - X \beta|_2^2 + \lambda\|D \beta\|_1,
#
# where D is (p-1) x p is the discrete difference operator (of
# order 1), and X is n x p with full column rank. The solution is
# piecewise constant, with adaptively chosen break points.

fusedlasso1d <- function(y, pos, X, gamma=0, approx=FALSE, maxsteps=2000,
                         minlam=0, tol=1e-11, verbose=FALSE) {
  if (missing(pos)) pos = NULL
  if (missing(X)) X = NULL
  if (!is.null(X) && !is.matrix(X)) stop("X must be a matrix.")
  if (!is.null(pos) && !is.null(X)) stop("Underlying positions cannot be specified with a (nonidentity) predictor matrix X.") 
  if (is.null(X)) approx = TRUE
  
  if (!is.null(X) || gamma!=0) {
    nv = if(is.null(X)) length(y) else ncol(X)
    D = getD1dSparse(nv)
    out = fusedlasso(y,X,D,NULL,gamma,approx,maxsteps,minlam,tol,verbose,FALSE)
    out$ord = 0
    out$pos = pos
  }
  else out = trendfilter(y,pos,NULL,0,approx,maxsteps,minlam,tol,verbose)

  out$call = match.call()
  class(out) = c("fusedlasso", "trendfilter", "genlasso", "list")
  return(out)
}
