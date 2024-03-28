##-------------------------------------------------------------
##  PRplot绘制
##-------------------------------------------------------------
rm(list=ls())
library(PRROC)
library(SPARK)

itau <- 35
ipt <- 4
ifc <- 2
ipow <- 1
AN <- c(1, 5)#选择我们想要的策略
MethodName <- c("SPARK", "Mean", "0.05Add+0.1Per", "0.05Add+0.5minPer", "SPARK-Addition", "0.5minAdd+0.5minPer")
labels <- c(rep(1,1000), rep(0, 9000))
roc_data <- list()
ptn <- c("I", "II", "III")

for (cnt in 1:2) {
  an <- AN[cnt]
  if (an == 1) {
    dataname <- paste0("./sim_data/power_test/spark1/sim_MOB_power_spark1_pt", ipt,
                       "_fc", ifc, "_tau", itau, "_count")
    PRplotname <- paste0("./Plots/PRplot/PRplot_SPARK_pt", ipt, "_fc", ifc, "_tau", itau, "_count")
    ROCplotname <- paste0("./Plots/ROCplot/ROCplot_SPARK_pt", ipt, "_fc", ifc, "_tau", itau, "_count")
  } else {
    dataname <- paste0("./sim_data/power_test/spark2/sim_MOB_power_spark2_pt", ipt, 
                       "_fc", ifc, "_tau", itau, "_count")
    PRplotname <- paste0("./Plots/PRplot/PRplot_SPARK-Addition_pt", ipt, "_fc", ifc, "_tau", itau, "_count")
    ROCplotname <- paste0("./Plots/ROCplot/ROCplot_SPARK-Addition_pt", ipt, "_fc", ifc, "_tau", itau, "_count")
  }# 根据选用的不同分析方法设定不同的文件名/图片名
  load(paste0(dataname, ipow, ".RData"))
  
  auc_pr <- c()
  auc_roc <- c()
  
  
  # 根据我们选择的策略, 对pval进行不同的组合
  if (an > 2) {
    #  权重修改
    weights <- matrix(0,nrow=nrow(result),ncol=15)
    if (an == 3) {# SPARK3 加法核0.05, 周期核0.1
      for (i in 1:nrow(result)){
        for (j in 1:10)
          weights[i,j]=0.05
        for (j in 11:15)
          weights[i,j]=0.1
      }
    } else if (an == 4) {# SPARK4 加法核0.05, 周期核取最小
      for (i in 1:nrow(result)){
        for (j in 1:10)
          weights[i,j]=0.05
        weights[i, which.min(result[i,11:15])+10]=0.5
      }
    } else if (an == 5) {# SPARK5 所有核取最小
      for (i in 1:nrow(result)) {
        weights[i,which.min(result[i,1:15])]=1
      }
    } else if (an == 6) {# SPARK6 加法核取最小, 周期核取最小
      for (i in 1:nrow(result)){
        weights[i,which.min(result[i,1:10])]=0.5
        weights[i,which.min(result[i,11:15])+10]=0.5
      }
    } else if (an == 7) { # 取最小的四个
      for (i in 1:nrow(result)) {
        rnk <- rank(result[i,])
        for (j in 1:15) {
          if (rnk[j]<=4) {
            weights[i,j] <- 0.25
          }
        }
      }
    }
    res_pval=result[,1:15]
    combined_pvalue <- CombinePValues(res_pval, weights=weights)
    result <- data.frame(res_pval, combined_pvalue = combined_pvalue,  adjusted_pvalue = p.adjust(combined_pvalue, method="BY") )
  }#  权重修改
  
  # 以下都是对pval数据的处理, 以化为pr\roc函数能接收的格式
  if (an == 1) {
    pval <- cbind(result[,12])
  } else {
    pval <- cbind(result[,17])
  }
  rownames(pval) <- rownames(result)
  gene <- rownames(pval)
  for (row in 1:10000) {
    gene[row] <- gsub('[gene]', '', gene[row])
  }
  gene <- as.numeric(gene)
  pval <- cbind(pval, gene)
  pval <- cbind(pval[order(pval[,2]),1])
  pval[,1] <- 1 - pval[,1]
  
  # 下面开始pr\roc曲线的绘制以及曲线下面积的计算
  # pr数据的计算
  roc_data <- append(roc_data, list(roc.curve(scores.class0 = pval, weights.class0 = labels, curve=T)))
}

# pr曲线的绘制

png(paste0(ROCplotname, ipow, "(1).png"), width = 8*600, height = 8*600, res = 600)

plot(roc_data[[1]], col = "blue", xlab = "FDR", ylab = "Sensitivity", lty = 1, lwd = 5, 
     main = paste0("ROC - Curve", " (Pattern ", ptn[ipt-1], ")"), cex.main = 1.5, auc.main = FALSE)
plot(roc_data[[2]], col = "red", lty = 1, lwd = 5, add = TRUE)
legend("bottomleft", legend = c(MethodName[1], MethodName[5]), col = c("blue", "red"), lty = c(1, 1), lwd = 2)

text(0.45, 0.03, paste(MethodName[1], " AUC:", round(roc_data[[1]]$auc, 7)), col = "blue")
text(0.5, 0, paste(MethodName[5], " AUC:", round(roc_data[[2]]$auc, 7)), col = "red")

dev.off()



# pr曲线的保存
# png(paste0(PRplotname, ipow, ".png"), width = 8*600, height = 8*600, res = 600) # 将 res 参数设置为 300，表示图像的精度为 300 dpi
# svg(paste0(PRplotname, ipow, ".svg"))
# plot(pr_data, main = paste0("PR Curve - ", MethodName[an]))
# dev.off()
# p <- plot(pr_data, main = paste0("PR Curve - ", MethodName[an]))
# ggsave(filename = paste0(PRplotname, ipow, ".png"), plot = p, width = 2000, height = 2000, dpi = 300)



