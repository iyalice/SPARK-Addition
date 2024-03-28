
##-------------------------------------------------------------
## Spatial Distribution of Representative Genes
## Please set the work directory to "C:\\Users\\bobo\\Desktop\\DaChuang_Files\\My_work\\02-Data_Visualization"
##-------------------------------------------------------------

rm(list = ls())
library("readxl")
source("./funcs/funcs.R")

target = 'Rep11_MOB'
file_name = paste0("./Plots/Feature_plot/FeaturePlot_", target, ".png")

#target可在
#Rep11_MOB 小鼠嗅球切片11
#Rep12_MOB 小鼠嗅球切片12
#Layer2_BC 人体乳腺癌切片2
#之间切换

#导入基因读数
counts <- read.table(paste("./raw_data/", target, "_count_matrix-1.tsv", sep = "", collapse = NULL), check.names = F)

# 找到匹配的行索引
rows_to_extract <- which(row.names(result) %in% result21)
# 通过索引提取相应的行
extracted_rows <- result[rows_to_extract, ]
# 根据 column_to_sort 列的值从高到低获取排序后的行索引顺序
sorted_index <- order(extracted_rows$adjusted_pvalue, decreasing = FALSE)
# 使用排序后的索引顺序重新排列数据框
sorted_data <- extracted_rows[sorted_index, ]
# 提取排名前十的行名
top_ten_row_names <- row.names(sorted_data)[1:10]
# 提取前十行的特定列数据
column_data <- sorted_data[1:10, 17]
# 将该列的数据格式化为科学计数法（保留3位小数）
formatted_column_data <- format(column_data, scientific = TRUE, digits = 3)
# 使用 paste0() 合并行名与科学计数法的表达式，并用空格隔开
merged_data <- paste0(top_ten_row_names, " ", formatted_column_data)


#输入将要绘制的基因的名称, 貌似不能只画一个基因
gene <- top_ten_row_names
gene_plot <- gene

rn <- rownames(counts)

info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), "[", 1)),
                         y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

vst_ct <- var_stabilize(t(counts)) # R function in funcs.R
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)

# 每一个feature plot的title
genetitle <- merged_data

pltdat <- cbind.data.frame(info[,1:2],rel_vst_ct)

pp <- lapply(1:(ncol(pltdat)-2),
             function(x){pattern_plot2(pltdat,x,main=T,titlesize=1.5,title=genetitle[x])})
grid.arrange(grobs = pp, ncol = 5)

# 绘制
g <- arrangeGrob(grobs = pp, ncol = 5)

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






