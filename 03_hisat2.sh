#!/bin/bash

# 1. 定义变量
INDEX_PREFIX="./ref/index/IRGSP-1.0" # 索引文件路径
CLEAN_DATA="./results/clean_data"  # 注意：这是01_fastp.sh的输出目录
OUT_DIR="./results/sam"  # 输出路径
THREADS=4

# 2. 创建输出目录
mkdir -p ${OUT_DIR}

# 3. 循环读取 samples.txt 进行比对
for sample in $(cat samples.txt)
do
    echo "正在比对样本: ${sample} ..."

    hisat2 -x ${INDEX_PREFIX} \
           -1 ${CLEAN_DATA}/${sample}_1_clean.fastq.gz \
           -2 ${CLEAN_DATA}/${sample}_2_clean.fastq.gz \
           -S ${OUT_DIR}/${sample}.sam \
           -p ${THREADS}
done

echo "所有样本比对完成！"