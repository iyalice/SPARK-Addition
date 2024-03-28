
##-------------------------------------------------------------
## Spatial Distribution of Representative Genes
##-------------------------------------------------------------

rm(list = ls())
library("readxl")
source("./funcs/funcs.R")

ipt <- 2
ifc <- 2
itau <- 20
target = paste0("pt", ipt, "_fc", ifc, "_tau", itau)
file_name = paste0("./Plots/Feature_plot/SimData/pattern", ipt, "/FeaturePlot_", target, ".png")


#导入基因读数
counts <- read.table(paste("./sim_data/power_test/data/pattern", ipt,"/sim_MOB_power_", target, "_count1.csv", sep = "", collapse = NULL), check.names = F)

# 数据处理
data_split <- data.frame(do.call(rbind, strsplit(as.character(counts$V2), ",")))
counts <- cbind(counts, data_split)
counts <- counts[, -c(2,3)]
rownames(counts) <- counts[, 1]
counts <- counts[, -1]
colnames(counts) <- unlist(counts[1, ])
counts <- counts[-1,]

#输入将要绘制的基因的名称, 貌似不能只画一个基因
gene <- "gene1"
gene_plot <- gene

rn <- rownames(counts)

info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), "[", 1)),
                         y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

vst_ct <- var_stabilize(t(counts)) # R function in funcs.R
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)

# 每一个feature plot的title
# genetitle <- merged_data

pltdat <- cbind.data.frame(info[,1:2],rel_vst_ct)

# pp <- lapply(1:(ncol(pltdat)-2),
#              function(x){pattern_plot2(pltdat,x,main=T,titlesize=1.5,title=genetitle[x])})
pp <- lapply(1:(ncol(pltdat)-2),
             function(x){pattern_plot2(pltdat,x,main=T)


grid.arrange(grobs = pp, ncol = 1)

# 绘制
g <- arrangeGrob(grobs = pp, ncol = 1)

# 保存
ggsave(file = paste('./Plots/Feature_plot/Best_performance_gene/', target, '.png', sep = "", collapse = NULL), width = 297, height = 210, units = "mm", g)

#---------------------------------------------------------
# 以下是批量绘制图片

# 将数据转为excel
# mydata <- read_excel(paste('./Result/Result_of_', target, '.xlsx', sep = "", collapse = NULL))

# for (i in 1:nrow(mydata)) {
#   
#   #数据预处理
#   gene <- mydata[i,]
#   if (i == nrow(mydata)) { # 最后的数据可能有空值
#     gene <- t(na.omit(t(gene)))
#   }
#   gene <- gsub('["]', '', gene) # 去掉原数据中的双引号
#   gene <- gsub("\r", '', gene) # 去掉原数据中可能出现的"\r"
#   gene <- gsub("\n", '', gene) # 去掉原数据中可能出现的"\n"
#   
#   #绘图
#   gene_plot <- c(gene)
#   
#   rn <- rownames(counts)
#   info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), "[", 1)),
#                            y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
#   rownames(info) <- rn
#   
#   vst_ct <- var_stabilize(t(counts)) # R function in funcs.R
#   sig_vst_ct <- vst_ct[gene_plot, ]
#   rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)
#   
#   genetitle <- c(gene)
#   pltdat <- cbind.data.frame(info[,1:2],rel_vst_ct)
#   
#   pp <- lapply(1:(ncol(pltdat)-2),
#                function(x){pattern_plot2(pltdat,x,main=T,titlesize=1.5,title=genetitle[x])})
#   grid.arrange(grobs=pp, ncol=5)
#   
#   g <- arrangeGrob(grobs = pp, ncol = 5)
#   
#   ggsave(file = paste(file_name, '[', i, '].png', sep = "", collapse = NULL), width = 297, height = 210, units = "mm", g)
# 
# }

#gene_plot <- c("Aldh7a1", "Aim1", "Dgkb", "Wisp1", "Bpifb9a", "Scd2", "Ttc3", "Lin7a", "Cav1", "Pja1")
#vst_ct <- var_stabilize(t(counts)) # R function in funcs.R
#sig_vst_ct <- vst_ct[gene_plot, ]
#rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)

#genetitle <- c(expression("Aldh7a1"))

