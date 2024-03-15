rm(list=ls())
library(SPARK)
library("data.table")
Rep12_MOB_count_matrix <- fread("~/Rep12_MOB_count_matrix-1.tsv")
col1<-Rep12_MOB_count_matrix$V1
as.character(col1)
row.names(Rep12_MOB_count_matrix) <- col1
Rep12_MOB_count_matrix <- Rep12_MOB_count_matrix[,-1]
Rep12_MOB_count_matrix <- t(Rep12_MOB_count_matrix)
colnames(Rep12_MOB_count_matrix)<- col1
info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(Rep12_MOB_count_matrix),split="x"),"[",1)),
                         y=as.numeric(sapply(strsplit(colnames(Rep12_MOB_count_matrix),split="x"),"[",2)),
                         total_counts=apply(Rep12_MOB_count_matrix,2,sum))
rownames(info) <- colnames(Rep12_MOB_count_matrix)
spark <- CreateSPARKObject(counts=Rep12_MOB_count_matrix, 
                           location=info[,1:2],
                           percentage = 0.1, 
                           min_total_counts = 10)
spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, 
                  covariates = NULL, 
                  lib_size = spark@lib_size, 
                  num_core = 5,
                  verbose = F)
X=spark@location[,1:2]
kmt=list()
ED <- as.matrix(dist(X))
lrang <- ComputeGaussianPL(ED, compute_distance=FALSE)[3:7]
a=1
for(ikernel in 2:5){
  for (jkernel in 1:(ikernel-1))
  {
    kernel_mat <- exp(-ED^2/(2*lrang[ikernel]^2))+exp(-ED^2/(2*lrang[jkernel]^2))
    kmt[[a]]<-kernel_mat
    rm(kernel_mat)
    a = a+1
  }
} 
for(ikernel in 1:5){
  kmt[[a]] <- cos(2*pi*ED/lrang[ikernel])
  a=a+1 
}
spark <- spark.test(spark, 
                    kernel_mat=kmt,
                    check_positive = T, 
                    verbose = F)
result=spark@res_mtest[order(spark@res_mtest$combined_pvalue,decreasing = F),]
res_pval=result[,1:15]
combined_pvalue <- apply(res_pval, 1, min)
result <- data.frame(res_pval, combined_pvalue = combined_pvalue,  adjusted_pvalue = p.adjust(combined_pvalue, method="BY") )
save(result,file = "~/Rep12_MOB.RData")
