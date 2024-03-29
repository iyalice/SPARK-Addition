# SPARK-Addition
SPARK-Addition是一个用于鉴别空间可变基因的程序，这是生物信息学的热门话题。我们在SPARK的基础上改进了呈现空间表达模式的核函数，并修改了P值的计算规则。这使得 SPARK-Addition 的性能比原始 SPARK 更好。

我们将所有数据放在./data中，原始代码放在./analysis中，输出放在./output中，绘图和相关生成代码放在./Plots中，模拟数据生成代码放在./simulation中，并且我们还给出了程序的复杂版本./plus 。注意：需要将gpml包置于./plus下，才能使main.R或optimizer_2.m正常运行。

# 程序使用示例
### 数据预处理与SPARK的构建
我们在示例中所使用的真实数据集为Layer2_BC_Count.rds，该数据为关于人乳腺癌空间转录组的基因表达计数矩阵。其中行名为基因名称，列名为二维坐标对信息，该矩阵中所有表达量均为整数。在预处理步骤中，SPARK-Addition需要提取二维坐标的横纵坐标值，并对每个坐标点上全部的基因表达求和，在本例中结果最终存储到名为info的数据框。“CreateSPARKObject”函数是生成spark的函数，我们使用了其中四个参数。现在对这四个参数做详细说明：counts接收的数据为p×n的基因表达计数矩阵，p是基因的数量而n是二维坐标点（或细胞）的数量；location接受的数据为表示位置信息（或是t-SNE/UMAP方法得到的两个成分）的n×2矩阵；percentage用于设置基因在全部二维坐标点中有表达量的占比的阈值；min_total_counts用于设置二维坐标点中所有的基因表达计数之和的阈值。本例选用Layer2_BC_Count.rds的基因表达计数矩阵rawcount（在读入该数据后即可得到）作为counts输入；info的前两列info[,1:2]为location输入；percentage分别设置为0.1和10。spark.vc拟合基于计数的空间模型来估计参数，除了必需提供的spark作为输入，我们同样调节了四个参数。Covariates提供了实验中的协变量，即混杂因素或批次效应的输入；lib_size用于获取每个二维坐标点的计数深度；num_core允许用户自定义拟合模型使用的内核数量；verbose则表示用户是否希望输出拟合信息。

![Image text](https://github.com/iyalice/SPARK-Addition/blob/main/code%20screenshot/Layer2_BC(1).png)

### 加法核的构建
X是spark中location变量所存储的二维坐标点拆分为横纵坐标值的信息，为n×2矩阵，ED对该矩阵中每个坐标对与其他坐标对的欧式距离进行计算，得到关于欧式距离的n×n实对称矩阵。ComputeGaussianPL函数将ED矩阵中所有的欧式距离进行对数变换，再按照数值的大小设置了10个区间分界点，而lrang存储了第3到第7个区间分界点。代码部分第一个for循环即为加法核的构建过程，也是SPARK-Addition的核心思路。即将lrang存储的5个点用于构建5个尺度参数不同的高斯核。加法核则在此之上将5个不同的高斯核进行两两加法组合，得到了10个全新的加法核，并记为kernel1-10。而第二个for循环用于构建5个尺度参数不同的周期核，记为kernel11-15。这15个核函数存储于kmt列表。

![Image text](https://github.com/iyalice/SPARK-Addition/blob/main/code%20screenshot/Layer2_BC(2).png)

### 核函数检验与结果筛选
spark.test函数根据提供的spark数据进行检验,并有三个可调节的参数。我们设置kernel_mat = kmt将存储在met列表中的10个加法核和5个周期核输入到函数之中；check_positive检测核函数矩阵是否为正；verbose则表示用户是否希望输出拟合信息。在经过spark.test函数检验后，需要对分析的数据进行筛选，这些数据储存于spark@res_mtest。该数据包括每个基因在每个核函数下的匹配程度，在数据中体现为大小不同的p值，每个基因在15个核中的最小p值记为combined_pvalue。对combined_pvalue使用Benjamini–Yekutieli程序控制错误发生率（FDR），因此得到人乳腺癌数据中每个基因的adjusted_pvalue。最终，我们筛选adjusted_pvalue<0.05的基因作为SPARK-Addition最终鉴定的SVG。

![Image text](https://github.com/iyalice/SPARK-Addition/blob/main/code%20screenshot/Layer2_BC(3).png)
