# 将表格导出成csv格式
# 清空内存
rm(list = ls())

# 加载保存的模块
load(file = 'table_data.Rdata')

# 1. 提取显著差异基因（padj < 0.05 且 变化幅度超过阈值）
# DEG_table[]:“[]”中表达式返回的是一个bool组成的向量，对应着每一个行index的筛选结果
# DEG_table[bool,]返回True的行
Data_choose <- DEG_table[DEG_table$padj < 0.05 & abs(DEG_table$log2FoldChange) > logFC_cutoff,]

# 2. 按 padj 从小到大排序（最显著的排最前面）
Data_choose <- Data_choose[order(Data_choose$padj),]

# 3. 把基因名字变成一列（因为 write.csv 默认不保留行名）
Data_choose$gene <- rownames(Data_choose)

# 4. 把gene变成第一列,setdiff(A,B)是一个做AB交集的函数，在A和B交集中的元素会被删除然后返回余集合按照原顺序
Data_choose <- Data_choose[, c("gene", setdiff(colnames(Data_choose), "gene"))]

# 5. 导出到 D 盘
# row.name参数,由于csv是二维表格，所以我们前边已经将行名变成第一列，这里不传入行名
write.csv(Data_choose,"D:/R/Rprogram/test/DEG_result.csv",row.names = FALSE)

# ========== KEGG/GO分析 ============
# 1. 安装 BiocManager（如果还没装的话）
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 2. 通过 BiocManager 安装 clusterProfiler 和必需的注释包
# GO（Gene Ontology）：基因本体论，主要讲基因的“功能”。
# KEGG（Kyoto Encyclopedia of Genes and Genomes）：京都基因与基因组百科全书，主要讲基因参与的“通路图”。
# clusterProfiler 是工具/引擎，它内部写了不同的“接口（API）”，可以分别去连接 GO 的服务器和 KEGG 的服务器。
BiocManager::install("clusterProfiler")

# airway 是人类数据，所以需要 org.Hs.eg.db（人类基因注释数据库）
# org.Hs.eg.db:是一个物种注释数据库（Annotation Database）。
#   org：Organism（物种）。
#   Hs：Homo sapiens（人类）。
#   eg：Entrez Gene ID（NCBI 的基因编号）。
#   .db：Database（数据库）。
# 为什么需要它？
# 因为 R 不知道人类基因长什么样。你需要给它一本“人类基因字典”，它才能把 ENSG00000152583 翻译成 SPARCL1，并告诉 R 这个基因对应哪个 GO 条目、哪个 KEGG 通路。
BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)

# 从 Data_choose 中提取基因 ID（Ensembl ID，如 "ENSG00000152583"）
# rownames(Data_choose) 返回的是带有版本号的 Ensembl ID
gene_list <- rownames(Data_choose)

# 或者，如果你只想分析上调基因：
# up_genes <- rownames(Data_choose[Data_choose$log2FC == "UP", ])

# 对全部差异基因进行 GO 富集分析
ego <- enrichGO(
  gene          = gene_list,     # 你的差异基因 ID 列表
  OrgDb         = org.Hs.eg.db,  # 基因对应的物种数据库
  keyType       = "ENSEMBL",     # 我的导入的矩阵中的gene的名称，如ENSG00000152583，是ENSEMBL格式 
  ont           = "ALL",         # 从三个维度进行基因分析 —— BP(生物学过程), CC(细胞组分), MF(分子功能)，得出结果时包括此gene在细胞什么位置，干什么活，参与了什么过程
  pAdjustMethod = "BH",          # BH = Benjamini-Hochberg，是一种多重检验校正方法。这个校正的不是传入表中的pvalue,而是GO分析结果的p值进行矫正
  pvalueCutoff  = 0.05,          # p 值阈值，并非传入矩阵的 p值阈值，而是GO分析后的p值阈值
  qvalueCutoff  = 0.2,           # q 值阈值，GO分析后的padj(BH多重检验校正的p值)阈值
  readable      = TRUE           # 将GO富集后return的Entrez ID转换成人类能看懂的gene ID
)

# 查看结果的前几行
head(ego)

ego_df <- data.frame(ego)


# 加载包
library(clusterProfiler)
library(org.Hs.eg.db)

# 1. 把 ENSEMBL ID 转换为 Entrez ID
entrez_ids <- clusterProfiler::bitr(
  gene_list, 
  fromType = "ENSEMBL", 
  toType   = "ENTREZID", 
  OrgDb    = org.Hs.eg.db
)


# 看一眼转换结果（转换过程中会丢掉一些没有对应 Entrez ID 的基因，这是正常的）
head(entrez_ids)

# 对全部差异基因进行 KEGG 通路富集分析
ekegg <- enrichKEGG(
  gene          = entrez_ids$ENTREZID,      # 基因 ID 列表，KEGG 默认接受 Entrez ID,Entrez ID是NCBI给每一个gene的唯一数字编号
  organism      = "hsa",          # 物种代码，"hsa" 代表人类
  keyType       = "kegg",         # 基因 ID 类型
  pAdjustMethod = "BH",           # 多重检验校正方法
  pvalueCutoff  = 0.05,           # p 值阈值
  qvalueCutoff  = 0.2             # q 值阈值
)

# 查看结果
head(ekegg)

kegg_df <- data.frame(ekeg)

library(ggplot2)

# 1. 绘制 GO 富集条形图（显示最显著的 20 个 GO 条目）
p1 <- barplot(ego, showCategory = 20, title = "GO Enrichment Analysis", label_format = NULL)
ggsave("D:/R/Rprogram/test/GO_barplot.png", plot = p1, width = 15, height = 10)

# 2. 绘制 GO 富集气泡图（更美观，显示 p 值和基因数量）
p2 <- dotplot(ego, showCategory = 20, title = "GO Enrichment Analysis")
ggsave("D:/R/Rprogram/test/GO_barplot2.png", plot = p2, width = 20, height = 10)

# 3. 绘制 KEGG 富集气泡图
p3 <- dotplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment")
ggsave("D:/R/Rprogram/test/KEGG_barplot.png", plot = p3, width = 20, height = 10)
