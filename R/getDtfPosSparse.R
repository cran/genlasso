getDtfPosSparse <- function(n,k,z) {
  weights = 1/diff(c(z,z[n]+1))
#  weights = weights/max(weights)
  D = bandSparse(n, m=n, c(0,1), diagonals=list(-weights,weights))
  D0 = bandSparse(n, m=n, c(0,1), diagonals=list(rep(-1,n),rep(1,n-1)))
#  D0 = bandSparse(n, m=n, c(0,1), diagonals=list(-weights,weights))
  if (k != 0) {
    for (i in 1:k) {
      D = D0 %*% D
#      D = D/max(D)
    }
  }
  return(D[1:(n-k-1),])
}
