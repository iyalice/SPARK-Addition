
##-------------------------------------------------------------
## Spatial Distribution of Representative Genes
## Please set the work directory to "C:\\Users\\bobo\\Desktop\\DaChuang_Files\\My_work\\02-Data_Visualization"
##-------------------------------------------------------------

rm(list = ls())
library("readxl")
source("./funcs/funcs.R")

target = 'Layer2_BC'
file_name = paste0("./Plots/Feature_plot/FeaturePlot_", target, "_MICe.png")

#target可在
#Rep11_MOB 小鼠嗅球切片11
#Rep12_MOB 小鼠嗅球切片12
#Layer2_BC 人体乳腺癌切片2
#之间切换

#导入数据
counts <- read.table(paste("./raw_data/", target, "_count_matrix-1.tsv", sep = "", collapse = NULL), check.names = F)


#输入将要绘制的基因的名称, 貌似不能只画一个基因
gene <- c("DUT", "ARHGAP23", "GLRX5", "EVA1B", "CALM3")
gene_plot <- gene

rn <- rownames(counts)

info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), "[", 1)),
                         y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

vst_ct <- var_stabilize(t(counts)) # R function in funcs.R
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)

# 每一个feature plot的title
genetitle <- gene

pltdat <- cbind.data.frame(info[,1:2],rel_vst_ct)

pp <- lapply(1:(ncol(pltdat)-2),
             function(x){pattern_plot2(pltdat,x,main=T,titlesize=1.5,title=genetitle[x])})
grid.arrange(grobs = pp, ncol = 5)

# 绘制
g <- arrangeGrob(grobs = pp, ncol = 5)

# 保存
ggsave(file = paste('./Plots/Feature_plot/Best_performance_gene/', target, '_MICe.png', sep = "", collapse = NULL), width = 297, height = 210, units = "mm", g)




