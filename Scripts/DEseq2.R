# # 环境配置和R包安装
# 
# # 1. 定义需要安装的包名（分为Bioconductor和CRAN两组）
# bioc_pkgs <- c("DESeq2", "apeglm", "clusterProfiler")
# cran_pkgs <- c("tidyverse", "pheatmap")
# 
# # 2. 检查Bioconductor包：如果其中任何一个没安装，则执行安装
# if (!all(sapply(bioc_pkgs, requireNamespace, quietly = TRUE))) {
#   BiocManager::install(bioc_pkgs)
# }
# 
# # 3. 检查CRAN包：如果其中任何一个没安装，则执行安装
# if (!all(sapply(cran_pkgs, requireNamespace, quietly = TRUE))) {
#   install.packages(cran_pkgs)
# }

library(DESeq2)

# c()代表 combine（合并）。
# 它的作用是把括号里用逗号隔开的几个元素，拼成一个“向量”（一串相同数据类型组成的列表）。
group_list=c('untrt','trt','untrt','trt','untrt','trt','untrt','trt') 

# colnames(Data):提取Data中的列名，也就是8个样本ID SRRXXX
# data.frame( ... ):是R中创建表格的函数;row.names和group_list是data.frame的两个参数
#   row.names = colnames(Data):将刚才提取出来的列名，变成表格的行名字
#   group_list = group_list:
colData <- data.frame(row.names = colnames(Data) , group_list = group_list)

# DESeqDataSetFromMatrix:R中的构造函数，作用是把零散的计数矩阵和分组信息打包成DESeq2能识别的对象
# countData参数:传入原始整数表达矩阵
# colData参数:传入样本分组表(行是样本,列是属性,这是一个分组标签，告诉DESeq2每一个样本的处理方式)
# design = ~ group_list:"~"代表"依赖于..."，~ group_list的意思是我要比较不同group_list标签之间的基因表达差异
# 因为group_list里有 trt 和 untrt，所以 DESeq2 就知道要去算这两组的差异基因。
dds <- DESeqDataSetFromMatrix(countData = Data,colData = colData,design = ~ group_list)

# 进行DESeq2的分析，经过三板斧的分析(归一化、估算离散度、拟合模型计算wald统计量)
#   离散度是同一个gene在不同样本间的离散度
ddsr <- DESeq(dds)

# results() 是 DESeq2 的一个函数,相当于一个static的函数,返回一个图表的对象。
# 它输入一个跑完 DESeq() 的 dds 对象，输出一张经过后处理（自动滤波、多重检验校正、可选差异倍数收缩）的统计结果表。
# contrast需要传一个向量参数，告诉results()“去 group_list 这一列里，用 "trt" 除以 "untrt"”
res <- results(ddsr,contrast = c("group_list","trt","untrt"))

# res1 <- as.data.frame(res)

# order() 函数不改变数值，它只返回一个行号索引（顺序）。由小到大。
# 它会找出 padj 列里最小的那个数在第几行，第二小的在第几行……然后返回c(3,1,2,4, ...) 这样的index number
# 然后res这个向量用return的索引来重新排列
resOrdered <- res[order(res$padj),]

head(resOrdered)

DEG_treat_2=as.data.frame(resOrdered)

# 去掉DEG中padj列中的NA行
DEG_treat_2 <- DEG_treat_2[!is.na(DEG_treat_2$padj),]

# copy这张表，作为画图的副本，保留原版的表格
DEG_table <- DEG_treat_2
