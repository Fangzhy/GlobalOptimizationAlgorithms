#Function of Thomson's problem, Fang Yu, 03/07/2016
Thomson.fn <- function(ThetaPhi) {
  fn.call <<- fn.call + 1
  ThetaPhi <- matrix(ThetaPhi, ncol = 2)
  xyz <- t(apply(ThetaPhi, 1, function(thetaphi) {
    c(sin(thetaphi[1]) * cos(thetaphi[2]),
      sin(thetaphi[1]) * sin(thetaphi[2]), cos(thetaphi[1]))}))
  #print(xyz)
  n <- nrow(ThetaPhi)
  tmp <- matrix(NA, nrow = n, ncol = n)
  index <- cbind(as.vector(row(tmp)), as.vector(col(tmp)))
  index <- index[index[, 1] < index[, 2], , drop=F]
  rdist <- apply(index, 1, function(idx) {
    tmp <- 1/sqrt(sum((xyz[idx[1], ] - xyz[idx[2], ])^2))})
  res <- sum(rdist)
  return(res)
}