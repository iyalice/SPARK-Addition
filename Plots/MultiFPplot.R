rm(list=ls())
library(SPARK)
library(ggplot2)
library(patchwork)
library(reshape2)
library(grid)
library(ggpubr)


# 选择要分析的数据
itau = 35 # 设定方差
ifc = 3 # 设定基因平均表达水平
N = 10000
maxFDR = 99
ptn <- c("I", "II", "III")
tau <- c(20, 35, 60)
# 所有方法的名字, 我们只选取第一和第五个
MethodName <- c("SPARK", "Mean", "0.05Add+0.1Per", "0.05Add+0.5minPer", "SPARK-Addition", "0.5minAdd+0.5minPer")
# plotname_png <- paste0("./Plots/FPplot_fc", ifc, "_tau", itau, "Allmin.png")
# plotname_svg <- paste0("./Plots/FPplot_fc", ifc, "_tau", itau, "Allmin.svg")
plot_list <- list()

# T/F: ord
# P/N: pvalue
TP <- matrix(0,nrow = maxFDR,ncol = 10)
FN <- matrix(0,nrow = maxFDR,ncol = 10)
power <- matrix(0,nrow = maxFDR,ncol = 10)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
vplayout = function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

# 循环, 绘制不同参数下的FDR曲线
for (itau in 1:3) {
  for (ifc in 2:4) {
    plotname_png <- paste0("./Plots/FDR_power_plot/FPplot_fc", ifc, "_tau", tau[itau], "Allmin.png")
    plotname_svg <- paste0("./Plots/FDR_power_plot/FPplot_fc", ifc, "_tau", tau[itau], "Allmin.svg")
    plot_list <- list()
    for (ipt in 2:4) {
      FDR <- c(1:maxFDR)
      for (i in 1:maxFDR) {FDR[i] = FDR[i]/100}
      FDR_power <- data.frame(FDR)
      for (an in 1:5) {
        
        if (an!=1 && an!=5) next
        
        analysis_by <- paste0("spark", an)
        if (an == 1) {
          dataname <- paste0("./sim_data/power_test/spark1/sim_MOB_power_spark1_pt", ipt, "_fc", ifc, "_tau", tau[itau], "_count")
        } else {
          dataname <- paste0("./sim_data/power_test/spark2/sim_MOB_power_spark2_pt", ipt, "_fc", ifc, "_tau", tau[itau], "_count")
        }
        for (count in 1:10){
          load(paste0(dataname, count, ".RData"))
          if (an > 2 && an < 7) {
            
            #  权重修改
            weights <- matrix(0,nrow=nrow(result),ncol=15)
            if (an == 3) {# SPARK3 加法核0.05, 周期核0.5
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
            }
            res_pval=result[,1:15]
            combined_pvalue <- CombinePValues(res_pval, weights=weights)
            result <- data.frame(res_pval, combined_pvalue = combined_pvalue,  adjusted_pvalue = p.adjust(combined_pvalue, method="BY") )
          } else if (an == 7) {
            res_pval=result[,1:15]
            min_pvalue <- minPValues(res_pval)
            result <- data.frame(res_pval, min_pvalue = min_pvalue,  adjusted_pvalue = p.adjust(min_pvalue, method="BH") )
          }
          
          for (fdr in 1:maxFDR) {
            result$ord <- as.numeric(gsub("[^0-9]","",rownames(result)))
            TP[fdr,count] = sum(result$ord<=1000 & result$adjusted_pvalue<=fdr/100)
            FN[fdr,count] = sum(result$ord<=1000 & result$adjusted_pvalue>fdr/100)
            power[fdr,count] = TP[fdr,count]/(TP[fdr,count]+FN[fdr,count])
          }
        }
        spark <- data.frame(rowMeans(power))
        colnames(spark) <- MethodName[an] # 设置曲线名字
        
        FDR_power <- cbind(FDR_power, spark)
      }
      FDR_power <- melt(FDR_power, id="FDR")
      colnames(FDR_power) <- c("FDR", "method", "Power")
      
      # create plot
      temp <- list(ggplot(data=FDR_power, aes(x=FDR, y=Power, group=method, color=method)) + 
                     geom_line(size = 2) + theme_bw() + coord_cartesian(ylim=c(0,1)) + 
                     labs(title = paste("Pattern", ptn[ipt-1])) +
                     theme(plot.title = element_text(family = "Lato", size=14, hjust=0.5),
                           legend.position = "none")
      )
      plot_list <- append(plot_list, temp)
    }
    
    # 最后一张图添加图例
    plot_list[[3]] <- plot_list[[3]] + theme(plot.title = element_text(family = "Lato", size=14, hjust=0.5),
                                             legend.position = c(0.7, .5),#plot内位置
                                             legend.text = element_text(family = "Lato", size = 10), # 设置图例文本大小
                                             legend.title = element_blank(),
                                             legend.key.size = unit(2, "lines"),
                                             # legend.title = element_text(family = "Lato", size = 12) # 设置图例标题大小)
    )
    # 将三个表达模式的图像整合到一起
    FPplot <- ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol = 3, nrow = 1)
    ggexport(filename = plotname_png, FPplot, width = 1500, height = 877, pointsize = 20)
    # ggexport(filename = plotname_svg, FPplot, width = 12.32, height = 5.6, pointsize = 20)
    print(FPplot)
   
  }
}

# 将所有FDR图像合并为一张
library(gridExtra)
library(png)
library(grid)
list_of_image_paths <- c()
for (itau in 1:3) {
  for (ifc in 2:4) {
    list_of_image_paths <- append(list_of_image_paths, paste0("./Plots/FDR_power_plot/FPplot_fc", ifc, "_tau", tau[itau], "Allmin.png"))
  }
}

imgs <- lapply(list_of_image_paths, function(path) {
  img <- readPNG(path)
  rasterGrob(img)
})

combined_img <- grid.arrange(grobs = imgs, ncol = 3)
plot(combined_img)


