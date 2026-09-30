#!/bin/bash
# ==============================================================================
# Pipeline: Salmon Index Construction & Quant Pipeline
# Author: Wuhaoyi
# Date: 2026-09-02
# Description: Standard Decoy-aware Salmon Quantification Pipeline
# ==============================================================================
 
set -e # 遇到错误立即停止执行
 
# 默认全局参考基因组与 CDS 路径
GENOME_FA="Ambrosia_artemisiifolia_genome.fa"
CDS_FA="Ambrosia_artemisiifolia_cds.fa"
THREADS=50
 
# 帮助信息函数
usage() {
    echo "Usage: $0 -i <fq_dir> -o <out_dir> -s <sample_list> [-p <threads>]"
    echo "  -i FASTQ 文件夹路径"
    echo "  -o 输出结果文件夹路径"
    echo "  -s 样本 ID 列表文件 (每行一个 SRR 号)"
    echo "  -p 线程数 (默认: 50)"
    exit 1
}
 
# 解析命令行参数
while getopts "i:o:s:p:h" opt; do
    case ${opt} in
        i) FQ_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        s) SAMPLE_LIST="$OPTARG" ;;
        p) THREADS="$OPTARG" ;;
        h|*) usage ;;
    esac
done
 
if [[ -z "$FQ_DIR" || -z "$OUT_DIR" || -z "$SAMPLE_LIST" ]]; then
    echo "[ERROR] 缺少必要参数！"
    usage
fi
 
# ==============================================================================
# 1. 检查或构建通用的 Decoy-aware Salmon Index
# ==============================================================================
INDEX_DIR="./global_salmon_index"
 
if [ ! -d "$INDEX_DIR" ]; then
    echo "========================================================================"
    echo "[STEP 1] 开始构建通用的 Salmon Decoy-aware Index..."
    echo "========================================================================"
    
    mkdir -p tmp_gentrome_dir
    
    # 官方标准步骤 1: 提取 Decoys (基因组染色体/Scaffold 标识符)
    grep "^>" "${GENOME_FA}" | cut -d " " -f 1 | sed 's/>//g' > tmp_gentrome_dir/decoys.txt
    
    # 官方标准步骤 2: 合并 CDS 与 Genome 生成 Gentrome
    cat "${CDS_FA}" "${GENOME_FA}" > tmp_gentrome_dir/gentrome.fa
    
    # 官方标准步骤 3: 建库 (-d 指定 decoys, --keepDuplicates 防止重复转录本被剔除)
    salmon index \
        -t tmp_gentrome_dir/gentrome.fa \
        -d tmp_gentrome_dir/decoys.txt \
        -i "${INDEX_DIR}" \
        -p "${THREADS}" \
        --keepDuplicates
        
    # 清理临时文件
    rm -rf tmp_gentrome_dir
    echo "[SUCCESS] Salmon Index 构建成功，存储于: ${INDEX_DIR}"
else
    echo "[INFO] 识别到已存在的 Salmon Index [${INDEX_DIR}]，跳过建库步骤。"
fi
 
# ==============================================================================
# 2. 批量定量样本 (Salmon Quant)
# ==============================================================================
mkdir -p "${OUT_DIR}"
 
echo "========================================================================"
echo "[STEP 2] 开始批量运行 Salmon 定量..."
echo "========================================================================"
 
while IFS= read -r sample || [[ -n "$sample" ]]; do
    # 过滤空行或以 # 开头的注释行
    [[ -z "$sample" || "$sample" =~ ^# ]] && continue
    
    # 自动跳过指定的离群样本
    if [[ "$sample" == "SRR10992850" || "$sample" == "SRR28099051" ]]; then
        echo "[SKIP] 自动跳过低比对率样本: ${sample}"
        continue
    fi
    
    # 匹配 FASTQ 文件路径 (兼容 .fq.gz 和 .fastq.gz)
    FQ1="${FQ_DIR}/${sample}_forward_paired.fq.gz"
    FQ2="${FQ_DIR}/${sample}_reverse_paired.fq.gz"
    
    if [ ! -f "$FQ1" ]; then
        # 尝试兼容部分中国数据的后缀命名差异 (.fq.g 等)
        FQ1=$(ls ${FQ_DIR}/${sample}*_1.f*q* ${FQ_DIR}/${sample}*forward*.f*q* 2>/dev/null | head -n 1)
        FQ2=$(ls ${FQ_DIR}/${sample}*_2.f*q* ${FQ_DIR}/${sample}*reverse*.f*q* 2>/dev/null | head -n 1)
    fi
 
    if [[ ! -f "$FQ1" || ! -f "$FQ2" ]]; then
        echo "[WARNING] 找不到样本 ${sample} 的 FASTQ 文件，跳过此样本！"
        continue
    fi
 
    echo ">>> 正在处理样本: ${sample} ..."
    
    salmon quant \
        -i "${INDEX_DIR}" \
        -l A \
        -1 "${FQ1}" \
        -2 "${FQ2}" \
        -o "${OUT_DIR}/${sample}_salmon_quant" \
        -p "${THREADS}" \
        --validateMappings \
        --gcBias \
        --seqBias
 
done < "${SAMPLE_LIST}"
 
echo "========================================================================"
echo "[SUCCESS] 所有样本定量完成！输出目录: ${OUT_DIR}"
echo "========================================================================"
 
 
