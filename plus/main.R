load("./plus/init.RData")
library("data.table")

gene_table = fread('./data/BC_pvals.csv')
gene_list = gene_table[,c("gene")]
X=spark@location[,1:2]
gene ="GAPDH"
cat(gene, file = './plus/current_gene.txt')
path=sprintf('./plus/%s.mat' ,gene)
path1=sprintf('./plus/%s_1.in' ,gene)
path2=sprintf('./plus/%s_2.in' ,gene)
path3=sprintf('./plus/%s.RData' ,gene)
y=spark@counts[c(gene),]
writeMat(path, X = X, y = y)
run_matlab_script('./plus/optimizer.m')
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
