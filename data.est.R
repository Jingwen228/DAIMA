
library(mvtnorm)
 test <- function(data1,interval=c(1,10),Time){
  
  phen <- data1$y2
  Time <- data1$Time
  geno <- data1$snp2
  
  nm <- dim(geno)[1]
  n1 <- interval[1]
  if(length(interval)==1){
    n2 <- interval[1]
  }else{ 
    n2 <- interval[2]}
  
  if(n2 >=nm)
  {n2 <- nm}
  res <- matrix(NA,nrow=length(c(n1:n2)),ncol=20)#存储LR,h02$value,以及8个参数
  for(i in n1:n2){
    SNP <- as.numeric(geno[i,])      #一个标记一个标记的求LR
    NSNP <- as.character(unlist(matrix(SNP,1)))
    missing <- which(NSNP==9)
    if ( length(missing) > 0)
    {
      SNP1 <- NSNP[-(missing)]#挑出第i个标记中未缺失的标记
      yy<- phen[-(missing), ]#然后挑出标记未缺失的表型
    }else{
      SNP1 <- NSNP
      yy <- phen
    }
    
    ndata <- data1
    ndata$y2 <- yy
    
    h01 <- try(lea.wei.H0 ( ndata),TRUE)#####通过调用后求得LH0
    if (class(h01) == "try-error") 
      h01 <- NA
    
 
    h02 <- try(lea.wei.H1(y11=yy,SNP1,init.par = h01$para,Time),TRUE)#####通过调用后求得LH1
    ######这样求得LH1-LH0才具有可比性。
    if (class(h02) == "try-error")
      h02 <- NA
    
    LR <- 2*(h01$value-h02$value)
    if(is.na(h01)||is.na(h02)){
      allpar <- c(LR,rep(NA,25))
    }else{
      allpar <- c(LR,h02$value,h02$para)
    }
    
    cat("snp", i, "=", allpar, "\n");
    res[(i-(n1-1)),(1:length(allpar))] <- allpar  #length(allpar)=10,即存储LR,h02$value,以及8个参数。
  }
  return(res)
  
}


lea.wei.H1 <- function(y11,SNP1,init.par,Time){
  
  
  lea.wei.Y <- y11
  
  index <- table(SNP1)
  Marker.type <- names(index)
  
  g.par <- c()
  Marker.index <- list()
  for(i in 1:length(Marker.type)){
    
    Marker.n <- which(SNP1==Marker.type[i])
    Marker.index[[i]] <- Marker.n
    g.par <- c(g.par,init.par[-(1:4)])
  }

  
  parinx <- c(init.par[1:4],g.par)
  rr <-  optim(parinx,fw.mle,lea.wei.Y=lea.wei.Y,Time=Time,lea.wei.X=Marker.index,ng=length(index),
               method="BFGS",control=list(maxit=10000))
  return(list(para=rr$par,value=rr$value))
}



fw.mle <- function(parinx,lea.wei.Y,lea.wei.X,Time,ng){
  
  len.cov <- 4
  len <- 0
  len.len <- 3
  par.covar <- parinx[1:len.cov]
  cov.mat <- AR1_get_matrix(par.covar,Time)
  par.curve <- parinx[-c(1:len.cov)]
  weight <- lea.wei.Y[,10]
  A1 <- c()
  for(i in 1:ng){
    mu <- c(LC.get.mu(par.curve[(1+len):(len+len.len)],Time[1:9]),mean(weight))
    ns <- length(lea.wei.X[[i]])
    y  <- lea.wei.Y[(lea.wei.X[[i]]),]
    Y.delt <-  y - matrix(rep(mu,ns),nrow=ns,byrow=T)
    
    pv <- dmvnorm( Y.delt, rep(0, NCOL(cov.mat)), cov.mat,log=T)
    
    A1 <- c(A1,-sum( pv ))
    len <- len +len.len
  }
  
  Lh1 <- sum(A1)
  Lh1
}
