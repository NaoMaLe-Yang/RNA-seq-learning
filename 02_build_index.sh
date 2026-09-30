#!/bin/bash

# 1、定义变量(相对路径)
REF_DIR = "./ref"
INDEX_DIR = "./ref/index"
THREADS = 4

# 解压下载的GTF格式文件
if test ! -f "${REF_DIR}/IRGSP-1.0_representative_transcript_exon_2026-02-05.gtf";
then
  echo "解压GTF文件"
  gunzip ${REF_DIR}/IRGSP-1.0_representative_transcript_exon_2026-02-05.gtf.gz
fi

# 提取剪接位点和外显子
echo "提取剪接位点和外显子..."
hisat2_extract_splice_sites.py ${REF_DIR}/IRGSP-1.0_representative_transcript_exon_2026-02-05.gtf > ${REF_DIR}/splicesites.tsv
hisat2_extract_exons.py ${REF_DIR}/IRGSP-1.0_representative_transcript_exon_2026-02-05.gtf > ${REF_DIR}/exons.tsv


# 构建 HISAT2 索引
echo "开始构建索引..."
hisat2-build -p ${THREADS} --ss ${REF_DIR}/splicesites.tsv --exon ${REF_DIR}/exons.tsv ${REF_DIR}/IRGSP-1.0_genome.fasta ${INDEX_DIR}/IRGSP-1.0

echo "索引构建完成！"