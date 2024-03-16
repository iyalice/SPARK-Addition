rm(list=ls())
source("./funcs/utilities.R")
library(SPARK)
library(patchwork)
library(ggplot2)
library(reshape2)
library(grid)
library(ggpubr)

# itau不等于35时result有10001行, 会出现下标出界的情况

itau = 60
perm_pvals <- c()
plot_list <- list()
N=10001
ptn <- c("I", "II", "III")
AN <- c(1, 5)
MethodName <- c("SPARK", "Mean", "0.05Add+0.1Per", "0.05Add+0.5minPer", "SPARK-Addition", "0.5minAdd+0.5minPer")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))
vplayout = function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
# 这里的文件位置均需要修改
dataname1 <- paste0("./sim_data/null_test/spark1/sim_MOB_null_spark1_tau", itau, "_count")
dataname2 <- paste0("./sim_data/null_test/spark2/sim_MOB_null_spark2_tau", itau, "_count")
plotame <- paste0("./Plots/QQplot/sim_MOB_null_tau", itau, "_count")

for (cnt in 1:2) {
  perm_pvals <- c()
  an <- AN[cnt]
  if (an != 1 && an != 5) { next }
  for (ipow in 1:10) {
    if (an==1) {
      load(paste0(dataname1, ipow, ".RData"))
      res_pval=result[,1:10]
      combined_pvalue <- CombinePValues(res_pval)
    } else {
      load(paste0(dataname2, ipow, ".RData"))
      if (an==2) {
        res_pval=result[,1:15]
        combined_pvalue <- CombinePValues(res_pval)
      } else if (an < 7) {
        weights <- matrix(0, nrow=N, ncol=15)
        if (an==3) {# SPARK3 加法核0.05, 周期核0.5
          for (i in 1:nrow(result)){
            min1=min(result[i,11:15])
            for (j in 1:10)
              weights[i,j]=0.05
            for (j in 11:15)
              weights[i,j]=0.1
          }
        } else if (an==4) {# SPARK4 加法核0.05, 周期核取最小
          for (i in 1:nrow(result)){
            for (j in 1:10)
              weights[i,j]=0.05
            weights[i, which.min(result[i,11:15])+10]=0.5
          }
        } else if (an==5) {# SPARK5 所有核取最小
          for (i in 1:nrow(result)) {
            weights[i,which.min(result[i,1:15])]=1
          }
        } else if (an==6) {# SPARK6 加法核取最小, 周期核取最小
          for (i in 1:nrow(result)){
            weights[i,which.min(result[i,1:10])]=0.5
            weights[i,which.min(result[i,11:15])+10]=0.5
          }
        }
        res_pval=result[,1:15]
        combined_pvalue <- CombinePValues(res_pval, weights=weights)
      } else if (an == 7) {
        res_pval=result[,1:15]
        combined_pvalue <- minPValues(res_pval)
        # result <- data.frame(res_pval, min_pvalue = min_pvalue,  adjusted_pvalue = p.adjust(min_pvalue, method="BH") )
      }
    }
    perm_pvals  <- cbind(perm_pvals, combined_pvalue)
  }
  
  #绘制QQplot
  temp <- list(qplot_gg(perm_pvals, pt.size=3, ax.txt.size=70, ax.title.size=70, ) + 
                 labs(title = MethodName[an]) + 
                 theme(plot.title = element_text(size =90, family = "serif", face = "bold", hjust = 0.5))
  )
  plot_list <- append(plot_list, temp)
  
}

# 绘图并保存
multiqqplot <- ggarrange(plotlist = plot_list, ncol = cnt, nrow = 1)
ggexport(filename = paste0("./Plots/QQplot/minQQplot_tau", itau, ".png"), multiqqplot, width = 4500, height = 2280, pointsize = 60)
# ggexport(filename = paste0("minQQplot_tau", itau, ".svg"), multiqqplot, width = 10, height = 5, pointsize = 20)
# plot(multiqqplot)







# 将所有FDR图像合并为一张
library(gridExtra)
library(png)
library(grid)
library(ggpubr)
list_of_image_paths <- c()
tau <- c(20, 35, 60)
for (itau in 1:3) {
  list_of_image_paths <- append(list_of_image_paths, paste0("./Plots/QQplot/minQQplot_tau", tau[itau], ".png"))
}

imgs <- lapply(list_of_image_paths, function(path) {
  img <- readPNG(path)
  rasterGrob(img)
})

combined_img <- grid.arrange(grobs = imgs, ncol = 1)
plot(combined_img)

# 保存需手动在右侧的图片展示区导出

# ggexport(filename = "./Plots/QQplot/AllQQplot.png", combined_img, width = 2000, height = 3000, pointsize = 20)


