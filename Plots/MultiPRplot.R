##-------------------------------------------------------------
##  PRplot绘制
##-------------------------------------------------------------
rm(list=ls())
library(PRROC)
library(SPARK)

AN <- c(1, 5)#选择我们想要的策略
MethodName <- c("SPARK", "Mean", "0.05Add+0.1Per", "0.05Add+0.5minPer", "SPARK-Addition", "0.5minAdd+0.5minPer")
labels <- c(rep(1,1000), rep(0, 9000))
pr_data <- list()

ptn <- c("I", "II", "III")
tau <- c(20, 35, 60)

aucpr_data <- array(data = NA, dim = c(2,3,3,3,10))


for (i1 in 1:2) { cnt <- i1 # 两种方法
  for (i2 in 1:3) { ifc <- i2+1 # 三种表达模式
    for (i3 in 1:3) { itau <- tau[i3] # 三种方差
      for (i4 in 1:3) { ifc <- i4+1 # 三种表达水平
        for (i5 in 1:10) { ipow <- i5 # 十次平行重复试验
          an <- AN[cnt]
          # 根据选用的不同分析方法设定不同的文件名/图片名
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
          } 
          
          load(paste0(dataname, ipow, ".RData"))
          
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
          } # 权重修改
          
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
          pr_data <- pr.curve(scores.class0 = pval, weights.class0 = labels, curve=T)
          aucpr_data[i1,i2,i3,i4,i5] <- pr_data$auc.integral
          
        }
      }
    }
  }
}


# ----------------
# 绘制箱线图
# ----------------
# 两种方法
# 三种表达模式
# 三种方差
# 三种表达水平
# 十次平行重复试验

# ipt <- 2
# itau <- 35
# ifc <- 2
# 
# for (ipt in 2:4) {
#   for (itau in 1:3) {
#     for (ifc in 2:4) {
#       plotname <- paste0("./Plots/boxplot/boxplot_pt", ipt, "_tau", itau, "_fc", ifc, ".png")
#       
#       # 绘制箱线图
#       data <- list(
#         SPARK <- aucpr_data[1,ipt-1,itau,ifc-1,],
#         SPARK_Addition <- aucpr_data[2,ipt-1,itau,ifc-1,]
#       )
#       names(data) <- c("SAPRK", "SPARK-Addition")
#       
#       png(plotname, width = 8*600, height = 8*600, res = 600)
#       
#       boxplot(data,
#               xlab = "Method", ylab = "AUC_PR",
#               main = paste0("Pattern ", ptn[ipt-1], " (tau = 0.", tau[itau], " fc = ", ifc, ")"),
#               col = c("red", "blue"))  # 设置箱线图的颜色
#       
#       dev.off()
#     }
#   }
# }
for (ipt in 2:4) {
  for (ifc in 2:4) {
    plotname <- paste0("./Plots/boxplot/boxplot_pt", ipt, "_fc", ifc, ".png")
    data <- list(
      SPARK_tau20 <- aucpr_data[1,ipt-1,1,ifc-1,],
      SPARK_Addition_tau20 <- aucpr_data[2,ipt-1,1,ifc-1,],
      SPARK_tau35 <- aucpr_data[1,ipt-1,2,ifc-1,],
      SPARK_Addition_tau35 <- aucpr_data[2,ipt-1,2,ifc-1,],
      SPARK_tau60 <- aucpr_data[1,ipt-1,3,ifc-1,],
      SPARK_Addition_tau60 <- aucpr_data[2,ipt-1,3,ifc-1,]
    )
    names(data) <- c("SAPRK_tau20", "SPARK-Addition_tau20", "SAPRK_tau35", "SPARK-Addition_tau35", "SAPRK_tau60", "SPARK-Addition_tau60")
    names(data) <- c("SAPRK", "SPARK-Addition", "SAPRK", "SPARK-Addition", "SAPRK", "SPARK-Addition")
    png(plotname, width = 8*600, height = 8*600, res = 600)
    
    boxplot(data,
            xlab = "TAU = 0.20                         TAU = 0.35                         TAU = 0.60",
            ylab = "AUC_PR",
            main = paste0("Pattern ", ptn[ipt-1], " (fc = ", ifc, ")"),
            col = c("blue", "red"))  # 设置箱线图的颜色
    
    dev.off()
  }
}

library(gridExtra)
library(png)
library(grid)
list_of_image_paths <- c()
for (ifc in 2:4) {
  for (ipt in 2:4) {
    list_of_image_paths <- append(list_of_image_paths, paste0("./Plots/boxplot/merge/boxplot_pt", ipt, "_fc", ifc, ".png"))
  }
}

imgs <- lapply(list_of_image_paths, function(path) {
  img <- readPNG(path)
  rasterGrob(img)
})

combined_img <- grid.arrange(grobs = imgs, ncol = 3)
plot(combined_img)







