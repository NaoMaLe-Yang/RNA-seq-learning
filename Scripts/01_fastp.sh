#!/bin/bash
# Fastq文件清洗(第一行是Hashbang，写这个封装代码运行所需解释器的绝对地址)

# 1. 定义变量（使用相对路径）
RAW_DATA="./data"
CLEAN_DATA="./results/clean_data"
QC_DIR="./results/fastqc_clean"
THREADS=4  # 核心数

# 2. 创建结果目录（强烈建议把 FastQC 的目录也创建了）
mkdir -p ${CLEAN_DATA}
mkdir -p ${QC_DIR}

# 3. 循环读取 samples.txt 里的样本名
for sample in $(cat samples.txt)
do
    echo "正在处理样本: ${sample} ..."

    fastp -i ${RAW_DATA}/${sample}_1.fastq.gz \
          -I ${RAW_DATA}/${sample}_2.fastq.gz \
          -o ${CLEAN_DATA}/${sample}_1_clean.fastq.gz \
          -O ${CLEAN_DATA}/${sample}_2_clean.fastq.gz \
          -h ${CLEAN_DATA}/${sample}.html \
          -j ${CLEAN_DATA}/${sample}.json \
          -w ${THREADS}

    # 顺手在循环里把 fastqc 也跑了
    fastqc ${CLEAN_DATA}/${sample}_1_clean.fastq.gz \
           ${CLEAN_DATA}/${sample}_2_clean.fastq.gz \
           -o ${QC_DIR} -t ${THREADS}
done

echo "所有样本 fastp 质控及 FastQC 检验完成！"