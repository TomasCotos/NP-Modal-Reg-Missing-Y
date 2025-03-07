require(nor1mix)


### Variable Aleatoria mixtura de dúas normais, a distancia entre das modas dependen dun valor x
#   Valores aleatorios
error.dep.x <- function(x, mu=c(-1,1), sig=c(0.25,0.25), w=c(0.5,0.5)) {
  n <- length(x)
  if (length(w)==1) w <- c(w,1-w)
  if (length(w)==2) w <- matrix(w,nrow=n,ncol=2,byrow=TRUE)
  er <- double(n)
  for (i in 1:n) {
    MW.20 <- norMix(mu=c(mu[1]+(mu[2]-mu[1])*x[i],mu[2] ),sigma=sig,w=w[i,], name="MW.20")
    er[i] <- rnorMix(1, MW.20)
  }
  return(er)
}
#   Densidade
den.error.dep.x <- function(x, v.x, mu=c(-1,1), sig=c(0.25,0.25), w=c(0.5,0.5)) {
  v <- data.frame(x=rep(x,rep(length(v.x),length(x))), v.x=rep(v.x,length(x)),den=NA)
  if (length(w)==1) w <- c(w,1-w)
  if (length(w)==2) w <- matrix(w,nrow=length(x),ncol=2,byrow=TRUE)

  for (i in 1:length(x)) {
    for (j in 1:length(v.x)) {
      MW <- norMix(mu=c(mu[1]+(mu[2]-mu[1])*x[i],mu[2]),sigma=sig,w=w[i,], name="MW")
      v$den[length(v.x)*(i-1)+j] <- dnorMix(v.x[j], MW)
    }
  }
  return(v)
}