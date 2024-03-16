##-------------------------------------------------------------
##  Analyze with SPARK in sim data
##-------------------------------------------------------------
rm(list = ls())

library(SPARK)

info <- read.csv(paste0("./data(1)/Rep11_MOB_info_spark.csv"), row.names = 1)

# 自己设置的核函数kmt
X=info[,1:2] #对于layer2BC行列数可能对不上
kmt=list()
ED <- as.matrix(dist(X))
lrang <- ComputeGaussianPL(ED, compute_distance=FALSE)[3:7]
a=1
for(ikernel in 2:5){ # 10个加法核
  for (jkernel in 1:(ikernel-1))
  {
    kernel_mat <- exp(-ED^2/(2*lrang[ikernel]^2))+exp(-ED^2/(2*lrang[jkernel]^2))
    kmt[[a]]<-kernel_mat
    rm(kernel_mat)
    a = a+1
  }
} 
for(ikernel in 1:5) # 5个周期核
{
  kmt[[a]]=cos(2*pi*ED/lrang[ikernel])
  a=a+1
}

for (analysis_by in 1:2) {
  for (ipow in 1:10) {
    countdata <- t(read.csv(paste0("./data(1)/sim_MOB_pattern2_fc3_tau35_count_null", ipow, ".csv"), row.names = 1))
    spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], 
                               percentage = 0.1, min_total_counts = 10)
    
    spark@lib_size <- info$total_counts
    
    # t1 <- proc.time()
    spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                      num_core = 1, verbose = T, fit.maxiter = 500)
    
    if (analysis_by == 1) {
      spark <- spark.test(spark, check_positive = T, verbose = T)  
    } else if (analysis_by == 2) {
      spark <- spark.test(spark, kernel_mat=kmt, check_positive = T, verbose = T) 
      # 跑15个核的时候需要加kmt
    }
    
    result=spark@res_mtest[order(spark@res_mtest$combined_pvalue,decreasing = F),]
    save(result,file=paste0("./data(1)/spark", analysis_by, "/sim_MOB_fc3_tau35_null", ipow, ".RData"))
    # time_comp <- proc.time() - t1
  }
}

