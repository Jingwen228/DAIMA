rm(list=ls())
library(mvtnorm)
source("data.lw.R")
source("data.curve.R")
source("data.est.R")


data1<- leaw.data1(genotype="m.snp.txt",
                   leaf ="leaf2.txt",
                   weight11="dry weight.txt" )



ret.H0 <- lea.wei.H0(data1)

ret.H1 <- test(data1,interval=c(11,15),Time)














