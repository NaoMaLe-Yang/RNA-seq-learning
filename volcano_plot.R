# 火山图R包
library(ggplot2)

# with()：临时进入 DEG_table 环境，列名可直接写，不用 DEG_table$log2FoldChange
# abs()：绝对值函数
# sd()：标准差函数
# 这行只是计算一个 logFC 阈值：所有基因 |log2FoldChange| 的均值 + 2 倍标准差
# 阈值算出来是 numeric，这里是 1.38198
# 后续筛选差异基因时，还要同时看 padj < 0.05：
#   上调：log2FoldChange >  logFC_cutoff
#   下调：log2FoldChange < -logFC_cutoff
logFC_cutoff <- with(DEG_table,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)))

# ifelse函数 : ifelse(条件, TRUE时从这个表达式的返回值中取对应行, FALSE时从这个表达式的返回值中取对应行)
# 首先，ifelse(函数1，函数2，函数3):ifelse会计算这三个函数，分别返回三个等长的向量(vector)，其中函数1只能返回逻辑向量(True or False)
# 然后，根据函数1返回的逻辑向量，True时从函数2返回的向量中取对应行的值，False时从函数3返回的向量中取对应行的值
DEG_table$log2FC <- as.factor(ifelse(
                                      DEG_table$padj < 0.05 & abs(DEG_table$log2FoldChange) > logFC_cutoff,
                                      ifelse(DEG_table$log2FoldChange > logFC_cutoff,'UP','DOWN'),
                                      'NOT')
                                      )

# nrow():number of row 总共存在多少行
# paste()：String拼接函数，\n为换行符
this_title <- paste('Cutoff for logFC is',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG_table[DEG_table$log2FC == 'UP',]),
                    '\nThe number of down gene is ',nrow(DEG_table[DEG_table$log2FC=='DOWN',])
                    )

# ggplot():新建画布,每一个"+"就是给画布新增一个配置
# 参数data:数据来源
# 参数aes:Aesthetic mapping（美学映射），把数据框里的列，映射到图上的视觉属性。
#   color:按标签列上色，区分 UP / DOWN / NOT
# geom_point()：告诉 ggplot “我要画散点图”。
#   alpha = 0.4：透明度。0 = 全透明，1 = 不透明。0.4 适合点很密集的时候，能看出哪里点叠得多，哪里稀疏。
#   size = 1.75：每个点的大小。
# ggtitle():把刚才用 paste() 拼好的字符串（包含阈值和上下调基因数），作为这幅图的大标题。
# theme(...)：修改图的具体细节。
#   plot.title = element_text(...)：专门针对标题的字体设置。
#   size = 15：标题字体大小（覆盖上面 base_size=20 的全局设置）。
#   hjust = 0.5：水平对齐方式，0 = 左对齐，0.5 = 居中，1 = 右对齐。这里是居中。
# scale_color_manual()：手动设置颜色映射。
#   values = c('blue','black','red')：给出三个颜色。
g <- ggplot(data=DEG_table,
            aes(x = log2FoldChange , y = -log10(pvalue) , color = log2FC)) + 
  geom_point(alpha = 0.4,size = 1.75) + 
  theme_bw(base_size = 20) + 
  ggtitle(this_title) + 
  theme(plot.title = element_text(size = 15,hjust = 0.5)) + 
  scale_color_manual(values = c('blue','black','red'))

print(g)

# 保存及其属性
ggsave(
  filename = "D:/R/Rprogram/test/volcano.png",
  plot = g,
  width = 8,
  height = 6,
  dpi = 300
)
