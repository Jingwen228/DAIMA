 

library(mvtnorm)
lea.wei.H0 <- function(data1){
  lea.wei.Y <-data1$y2 
  Time <- data1$Time
  
  parinx <- c( 0.5,1.917107,0.4,8.324484,11.55122,5.279945,0.5667437)
  rr <- optim(parinx,LW.mle,lea.wei.Y=lea.wei.Y,Time=Time,method="Nelder-Mead",control=list(maxit=10000))
  
  return(list(para=rr$par,value=rr$value))
  
}
LW.mle <- function(parinx,lea.wei.Y,Time){
  
  len.cov <- 4
  par.covar <- parinx[1:len.cov]
  par.curve <- parinx[-c(1:len.cov)]
  weight <- lea.wei.Y[,10]
  
  cov.mat <- AR1_get_matrix(par.covar,Time)
  mu <- c(LC.get.mu(par.curve,Time[1:9]),mean(weight))
  
  n <- dim(lea.wei.Y)[1]
  Y.delt <- lea.wei.Y-matrix(rep(mu,n),nrow=n,byrow=T)
  
  pv  <-dmvnorm( Y.delt,rep(0, NCOL(cov.mat)),cov.mat,log = T)
  # cat("pv",pv,"\n")
  Lh0 <- -sum(pv)
  Lh0
}
######用AR1求协方差矩阵###

AR1_get_matrix <- function(par, times, options=list()){
  m <- length(times)
  var_cov<-matrix(0,m,m)
  for (i in 1:m){
    for (j in 1:m){
      var_cov[i,j] <- par[1]^abs(j-i)*par[2]
      var_cov[i,10] <- par[3]^abs(10-i)*par[4]
      var_cov[10,j] <- par[3]^abs(10-j)*par[4]
    }
  }
  return(var_cov)}



LC.get.mu <- function(par, times, options=list())
{
  A <- par[1]/(1+par[2]*exp(-par[3]*times))
  return (A);
}
#$para
#[1]  0.9482906  2.7084707  0.1124620  9.1751007 11.2048398  4.6821383  0.5333978

#$value
#[1] 622.0977
