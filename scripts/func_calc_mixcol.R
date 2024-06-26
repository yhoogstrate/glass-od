#!/usr/bin/env R


mixcol_single <- function(c1, c2, ratio=0.5) {
  n <- names(c1)
  
  m1 <- col2rgb(c1)
  m2 <- col2rgb(c2)
  
  m <- (ratio * m2) + ((1 - ratio) * m1)
  r <- rgb(m[1,],m[2,],m[3,], maxColorValue = 255)
  names(r) <- n
  
  return(r)
}



mixcol <- function(c1, c2, ratio=c(0.5)) {
  out = c()
  
  for(rat in ratio) {
    out <- c(out, mixcol_single(c1, c2, rat))
  }
  
  return (out)
}


