covQTLs are a series of genes that coordinate and govern development of plant organs and whole biomass. Here, we expand a new mapping model based on functional mapping to detect covariation quantitative trait loci (QTLs), which, via a series of hypothesis tests, allows quantification of how QTLs regulate covariation between organ growth and biomass accumulation.

# load data and function
rm(list=ls())
library(mvtnorm)
source("data.lw.R")
source("data.curve.R")
source("data.est.R")


data1<- leaw.data1(genotype="m.snp.txt",
                   leaf ="leaf2.txt",
                   weight11="dry weight.txt" )


# null hypothesis and Alternative hypothesis

ret.H0 <- lea.wei.H0(data1)

ret.H1 <- test(data1,interval=c(1,10),Time)

