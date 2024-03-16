load("~/大创/init.RData")
library("data.table")
library(SPARK)
library(R.matlab)
library(matlabr)

gene_table = fread('~/大创/BC_pvals.csv')
gene_list = gene_table[,c("gene")]
X=spark@location[,1:2]
# (k in 1092:1152)
#{
  #gene = gene_list$gene[k]
  gene ="GAPDH"
  cat(gene, file = '~/大创/current_gene.txt')
  path=sprintf('~/大创/%s.mat' ,gene)
  path1=sprintf('~/大创/%s_1.in' ,gene)
  path2=sprintf('~/大创/%s_2.in' ,gene)
  path3=sprintf("~/大创/%s.RData" ,gene)
  
  y=spark@counts[c(gene),]
  writeMat(path, X = X, y = y)
  
  
  run_matlab_script('F:/Documents/大创/optimizer.m')
  KM1 <- matrix(,nrow=N,ncol=N)
  strker1 <- readLines(path1)
  KM2 <- matrix(,nrow=N,ncol=N)
  strker2 <- readLines(path2)
  for (i in 1:N)
    for (j in 1:i)
    {
      ed1=eval(parse(text=strker1))
      ed2=eval(parse(text=strker2))
      KM1[i,j]=KM1[j,i]=ed1
      KM2[i,j]=KM2[j,i]=ed2
    }
  
  spark <- spark.test(spark, 
                      kernel_mat=list(KM1,KM2),
                      check_positive = T, 
                      verbose = F)
  
  result=spark@res_mtest[order(spark@res_mtest$combined_pvalue,decreasing = F),]
  save(result,file = path3)
  gene_count = sum(result[,c("adjusted_pvalue")]<=0.05)
  cat(paste(gene,',',gene_count,'\n'),file = '~/大创/total_result.csv',append = TRUE)
  rm(KM1)
  rm(KM2)
  Sys.sleep(60)
#}
