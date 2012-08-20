getDtfSparse <- function(n,k) {
  D = bandSparse(n, m=n, c(0,1), diagonals=list(rep(-1,n),rep(1,n-1)))
  D0 = D
  if (k != 0) {
    for (i in Seq(1,k)) D = D0 %*% D
  }
  return(D[Seq(1,n-k-1),]) 
}
