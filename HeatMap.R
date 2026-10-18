# 热图
library(pheatmap)

choose_gene <- head(rownames(DEG_table),100)
# rownames(DEG_table)：从你的差异结果表中取出所有基因的 ID（比如 ENSG00000152583）。
# head(..., 100)：取前 100 个。因为你的 DEG_table 已经按 padj 排过序（从小到大），所以取出的就是 “最显著的 100 个差异基因”。

choose_matrix <- Data[choose_gene,]
# Data：你的原始表达矩阵,featureCounts结束后R直接读取的最原始的只包含reads counts的矩阵。
# [choose_gene, ]：按行名筛选，只把第 1 步选出的那 100 个基因的行拿出来。
# choose_matrix：此时你有了一个 100 行 × 8 列的小矩阵，里面是这 100 个基因在 8 个样本中的原始整数表达量。

choose_matrix=t(scale(t(choose_matrix)))
# t(choose_matrix)：转置矩阵，变成 8 行（样本） × 100 列（基因）。
# scale(...)：对每一列（每个基因）执行 Z-score 标准化
# 公式：(数值 - 该基因均值) / 该基因标准差

pheatmap(choose_matrix,filename = 'D:/R/Rprogram/test/DEG_top_100_heatmap.png')
# 上热图！！


