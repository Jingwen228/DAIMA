

leaw.data1 <- function(genotype="m.snp.txt",
                      leaf ="leaf2.txt",
                      weight11="dry weight.txt" ){
  
  r.snp <- read.table(genotype,header=T)
  r.leaf <- read.table(leaf,header=T)
  r.weight <- read.table(weight11,header=T)
  
  
  r.name <- colnames(r.snp)
  
  r.n1 <- c()
  for(i in 1:length(r.name)){
    r.n1 <- c(r.n1,as.numeric(strsplit(r.name[i],"X")[[1]][2]))
  } 
  
  s1 <- table(c(r.n1,r.leaf[,1],r.weight[,1]))
  label <- as.numeric(names(s1[which(s1==3)]))
  
  snp.s <- c()
  for(i in label){
    snp.s <- c(snp.s,which(r.n1==i))
  }
  
  leaf.s <- c()
  for(i in label){
    leaf.s <- c(leaf.s,which(r.leaf[,1]==i))
  }
  
  
  weight.s <- c()
  for(i in label){
    weight.s <- c(weight.s,which(r.weight[,1]==i))
  }
  
  
  snp2 <- r.snp[,snp.s]
  leaf2 <- r.leaf[leaf.s,-1]/10
  weight2 <- r.weight[weight.s,-1]*100
  
  y2 <- cbind(leaf2,weight2)
  
  Time <- c(1:10)
  list(snp2=snp2,y2=y2,Time=Time)
}


