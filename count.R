# 首先进行count.txt文件数据的读取，将其整理成矩阵

rm(list = ls()) # 清空变量

options(stringsAsFactors = F) # 让读入的数据类型都变成String类型

Data=read.table(
  "D:/Document/KeDa_Document/实验资料/Data_and_Image/Data analysis/Practice Project/Rice_genome/count/count.txt",
  header = TRUE, # 第一行列名不读入新矩阵
  )  # 从Linux的共享文件夹中读入数据

counts <- Data[,7:ncol(Data)] # 从原矩阵的第七列开始取一直到最后一列

rownames(counts) <- Data[,1] # 将原矩阵的第一列内容赋值给新矩阵当行index

