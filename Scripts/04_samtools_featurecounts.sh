#!/bin/bash

# =========================================================
# 04_samtools_featurecounts.sh
# 功能：SAM -> BAM -> 排序 -> 索引 -> flagstat
# 输入：03_hisat2.sh 输出的 results/sam/${sample}.sam
# 输出：results/bam/${sample}.sort.bam        (排序后 BAM)
#       results/bam/${sample}.sort.bam.bai    (索引)
#       results/flagstat/${sample}.flagstat.txt
# =========================================================

# 1. 定义变量
SAM_DIR="./results/sam"              # 03_hisat2.sh 的输出目录
BAM_DIR="./results/bam"              # 输出：BAM 和索引
FLAGSTAT_DIR="./results/flagstat"    # 输出：flagstat 统计
THREADS=4
SORT_MEM="768M"                      # 每个线程的内存上限，16G 机器保守设置

# 2. 创建输出目录
mkdir -p ${BAM_DIR}
mkdir -p ${FLAGSTAT_DIR}
mkdir -p ${OUT_DIR}

# 3. 循环读取 samples.txt 进行 SAM -> BAM 处理
for sample in $(cat samples.txt)
do
    echo "正在处理样本: ${sample} ..."

    # (1) SAM -> BAM（未排序，作为中间产物）
    samtools view -b -@ ${THREADS} \
        -o ${BAM_DIR}/${sample}.bam \
        ${SAM_DIR}/${sample}.sam

    # (2) 排序
    samtools sort -@ ${THREADS} -m ${SORT_MEM} \
        -o ${BAM_DIR}/${sample}.sort.bam \
        ${BAM_DIR}/${sample}.bam

    # (3) 建索引（默认生成 ${sample}.sort.bam.bai）
    samtools index -@ ${THREADS} ${BAM_DIR}/${sample}.sort.bam

    # (4) flagstat 质控统计(可以用MultiQC进行总体数据统计图表绘制)
    samtools flagstat -@ ${THREADS} ${BAM_DIR}/${sample}.sort.bam \
        > ${FLAGSTAT_DIR}/${sample}.flagstat.txt

    # (5) 删除中间产物，节省磁盘
    rm -f ${BAM_DIR}/${sample}.bam
done

echo "所有样本 SAM->BAM、排序、索引、flagstat 完成！"

# 4. 收集所有 BAM 文件路径（用数组，防止文件名带空格）
BAMS=()
for sample in $(cat samples.txt)
do
    BAMS+=("${BAM_DIR}/${sample}.sort.bam")
done

# 5. 一次性定量所有样本
echo "开始 featureCounts 定量，共 ${#BAMS[@]} 个样本..."
featureCounts \
    -T ${THREADS} \
    -a ${GTF} \
    -p --countReadPairs \
    -t exon -g gene_id \
    -o ${OUT_DIR}/gene_counts.txt \
    "${BAMS[@]}"

echo "featureCounts 定量完成！"
