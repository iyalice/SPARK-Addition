##------------------------------------------------------------
## 模拟数据的生成 Data Generation
## Generate the countdata based on the info from the realdata and
## patterns summarized from the SpatialDE result

## Fold Change
##-----------------
rm(list = ls())

# based_on <- "1" # spark1 or spark2

PorN <- "null" # power or null

# if (based_on == "1") { # 这一步load决定了我们用什么spark来生成模拟数据
#   load(paste0("./output/Rep11_MOB_spark.rds"))
# } else if (based_on == "2") {
#   load(paste0("./output/Rep11_MOB_spark2.rds"))
# }

load(paste0("./output/Rep11_MOB_spark.rds"))

dirname <- paste0("./sim_data/", PorN, "_test/data")

if (!dir.exists(dirname)){ # 生成文件夹
  dir.create(dirname) 
} else { 
  print(paste0(dirname, " already exists!")) 
}



info <- spark@location
info$total_counts <- spark@lib_size

beta <- sapply(1:length(spark@res_vc), function(x) {
  spark@res_vc[[x]]$coefficients
})
nb <- sapply(1:4, function(x) {
  log(x * exp(median(beta)))
})
tau2 <- sapply(1:length(spark@res_vc), function(x) {
  spark@res_vc[[x]]$theta[2]
})

load("./output/MOB_Pattern_SpatialDE.rds")

itau = 60 # 残差
if (PorN == "power") {
  numSignal <- 1000 # SVG的个数有1000个
} else if (PorN == "null") {
  numSignal <- 0 # SVG的个数有0个
}
numGene <- 10000 # 所有基因共有10000个
newN <- info$total_counts

for (ipt in 2:2) { # SE基因的三种表达模式
  
  dirname1 <- paste0(dirname, "/pattern", ipt)
  if (!dir.exists(dirname1)){ # 生成文件夹
    dir.create(dirname1) 
  } else { 
    print(paste0(dirname1, " already exists!")) 
  }
  for (ifc in 2:2) { # 乘数因子, 基因表达水平的比值
    
    pattern <- datlist[[paste0("pattern", ipt)]]
    grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
    uu <- c(nb[1], nb[ifc], nb[1])
    
    for (ipow in 1:10) { # 平行重复试验
      
      set.seed(ipow)
      beta0 <- uu[grp]
      lambda0 <- sapply(1:numSignal, function(x) {
        exp(beta0 + rnorm(length(beta0), 0, itau/100))
      })
      newCt0 <- lapply(1:numSignal, function(x) {
        rpois(length(lambda0[, x]), newN * lambda0[, x])
      })
      
      beta1 <- rep(uu[3], nrow(pattern))
      lambda1 <- sapply(1:(numGene - numSignal), function(x) {
        exp(beta1 + rnorm(length(beta1), 0, itau/100))
      }) 
      newCt1 <- lapply(1:(numGene - numSignal), function(x) {
        rpois(length(lambda1[, x]), newN * lambda1[, x])
      })
      
      countdata <- data.frame(rbind(do.call(rbind, newCt0), do.call(rbind, 
                                                                    newCt1)))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
      
      
      if (PorN == "Power") {
        write.csv(t(countdata), file = paste0(dirname1, "/sim_MOB_", PorN, "_pt", 
                                              ipt, "_fc", ifc, "_tau", itau, "_count", 
                                              ipow, ".csv"), 
                  row.names = T)
      } else {
        write.csv(t(countdata), file = paste0(dirname1, "/sim_MOB_", PorN, 
                                              "_tau", itau, "_count", 
                                              ipow, ".csv"),
                  row.names = T)
      }
    }
  }
}

write.csv(info, file = paste0(dirname, "/Rep11_MOB_info_spark.csv"), row.names = T)