rm(list=ls())
library(SPARK)
library(R.matlab)
library(matlabr)
load("./data/Layer2_BC_Count.rds")
rawcount[1:10,1:10]
info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",1)),
                         y=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",2)),
                         total_counts=apply(rawcount,2,sum))
rownames(info) <- colnames(rawcount)
spark <- CreateSPARKObject(counts=rawcount, 
                           location=info[,1:2],
                           percentage = 0.1, 
                           min_total_counts = 10)
spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, 
                  covariates = NULL, 
                  lib_size = spark@lib_size, 
                  num_core = 5,
                  verbose = F)
N=250


