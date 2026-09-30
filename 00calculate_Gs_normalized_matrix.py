import os
import sys
import argparse
import numpy as np
import pandas as pd
 
 
def parse_args():
    parser = argparse.ArgumentParser(description="Generate and save normalized sum matrix (Sum-then-Log strategy).")
    parser.add_argument("-t", "--tpm", required=True, help="Path to 00salmon_tpm_matrix.csv")
    parser.add_argument("-m", "--filtered_m", required=True, help="Path to filtered_tpm_matrix.csv")
    parser.add_argument("--psba-pos", required=True, help="Path to pos_psba_id.txt")
    parser.add_argument("--psba-all", required=True, help="Path to psba_id.txt")
    parser.add_argument("-l", "--lat", required=True, help="Path to sample_pop_Lat.csv")
    parser.add_argument(
        "-o", "--outdir", 
        default=" 00data/", 
        help="Output directory for saved matrix"
    )
    return parser.parse_args()
 
 
def calculate_geometric_mean(series):
    """计算非负序列的几何平均数 Gs"""
    pos_vals = series[series > 0]
    if len(pos_vals) == 0:
        return 0.0
    return np.exp(np.mean(np.log(pos_vals)))
 
 
def build_tx_to_gene_map(gene_ids):
    """构建基因与转录本的映射关系字典 (含补丁)"""
    tx_to_gene = {}
    for gid in gene_ids:
        if "evm.TU." in gid:
            tx_id = gid.replace("evm.TU.", "evm.model.")
        elif "novel_gene_" in gid:
            tx_id = gid.replace("novel_gene_", "novel_model_")
        else:
            tx_id = gid
        tx_to_gene[tx_id] = gid
 
    PATCH_DICT = {
        "evm.model.Hic_asm_15.1576.1.61aad7f8": "evm.TU.Hic_asm_15.1576",
        "evm.model.Hic_asm_2.2349.3.61aae6ca":  "evm.TU.Hic_asm_2.2349",
        "evm.model.Hic_asm_2.2351.1.61aae6cb":  "evm.TU.Hic_asm_2.2351",
        "evm.model.Hic_asm_8.427.1.61aaf263":   "evm.TU.Hic_asm_8.427"
    }
    for tx_patch, gid in PATCH_DICT.items():
        generic_tx = gid.replace("evm.TU.", "evm.model.")
        if generic_tx in tx_to_gene:
            del tx_to_gene[generic_tx]
        tx_to_gene[tx_patch] = gid
 
    return tx_to_gene
 
 
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
 
    # 1. 读取矩阵与元数据
    full_tpm = pd.read_csv(args.tpm, index_col=0, sep=None, engine='python')
    filtered_tpm = pd.read_csv(args.filtered_m, index_col=0)
    meta_df = pd.read_csv(args.lat, sep=None, engine='python')
 
    pos_psba_ids = set(line.strip() for line in open(args.psba_pos) if line.strip())
    all_psba_ids = set(line.strip() for line in open(args.psba_all) if line.strip())
    non_pos_psba_ids = all_psba_ids - pos_psba_ids
 
    all_filtered_genes = list(filtered_tpm.index)
    
    # 2. 提取 3 类目标基因列表
    target_pos_psba = [g for g in pos_psba_ids if g in all_filtered_genes]
    target_non_pos_psba = [g for g in non_pos_psba_ids if g in all_filtered_genes]
 
    gene_info_path = os.path.join(os.path.dirname(args.filtered_m), "gene_info.csv")
    if os.path.exists(gene_info_path):
        gene_info = pd.read_csv(gene_info_path)
        nlr_genes_from_info = set(gene_info[gene_info['Group'] == 'NLR']['Gene_ID'])
        target_nlr = [g for g in all_filtered_genes if g in nlr_genes_from_info]
    else:
        target_nlr = [g for g in all_filtered_genes if g not in all_psba_ids and "psbA" not in g]
 
    # 3. 构建转录本映射规则
    tx_to_gene = build_tx_to_gene_map(all_filtered_genes)
    gene_to_tx = {v: k for k, v in tx_to_gene.items()}
 
    pos_psba_txs = [gene_to_tx[g] for g in target_pos_psba if gene_to_tx.get(g) in full_tpm.index]
    non_pos_psba_txs = [gene_to_tx[g] for g in target_non_pos_psba if gene_to_tx.get(g) in full_tpm.index]
    nlr_txs = [gene_to_tx[g] for g in target_nlr if gene_to_tx.get(g) in full_tpm.index]
 
    # 4. 计算 Top 500 HKG 的样本几何平均数 Gs
    gene_means = full_tpm.mean(axis=1)
    gene_sds = full_tpm.std(axis=1)
    cand_mask = gene_means >= 10.0
    cand_cvs = gene_sds[cand_mask] / gene_means[cand_mask]
    hkg_genes = cand_cvs.sort_values().head(500).index.tolist()
 
    hkg_tpm_matrix = full_tpm.loc[hkg_genes]
    gs_per_sample = hkg_tpm_matrix.apply(calculate_geometric_mean, axis=0)
 
    # 5. 执行【先求和 TPM，再除以 Gs，最后取 log2(ratio + 1)】的新策略
    samples = full_tpm.columns.tolist()
    summary_list = []
 
    for sample_id in samples:
        gs_val = gs_per_sample[sample_id]
        
        # 5.1 受正选择 psbA 物理 Sum 相对表达
        tpm_sum_pos = full_tpm.loc[pos_psba_txs, sample_id].sum()
        norm_sum_pos = np.log2((tpm_sum_pos / gs_val) + 1.0)
        
        # 5.2 未受正选择 psbA 物理 Sum 相对表达
        tpm_sum_non_pos = full_tpm.loc[non_pos_psba_txs, sample_id].sum()
        norm_sum_non_pos = np.log2((tpm_sum_non_pos / gs_val) + 1.0)
        
        # 5.3 48 个 NLR 物理 Sum 相对表达
        tpm_sum_nlr = full_tpm.loc[nlr_txs, sample_id].sum()
        norm_sum_nlr = np.log2((tpm_sum_nlr / gs_val) + 1.0)
        
        summary_list.append({
            "Sample": sample_id,
            "Sum_psbA_pos": norm_sum_pos,
            "Sum_psbA_non_pos": norm_sum_non_pos,
            "Sum_NLR": norm_sum_nlr
        })
 
    summary_df = pd.DataFrame(summary_list)
    final_df = pd.merge(summary_df, meta_df, on="Sample", how="inner")
 
    # 6. 保存矩阵文件
    out_file = os.path.join(args.outdir, "normalized_summary_matrix.csv")
    final_df.to_csv(out_file, index=False)
 
    # ================== Self-Check 检查与汇报 ==================
    print("========================================================================")
    print(f"[Self-Check] 成功更新标准化矩阵 (先求和后Log): {out_file}")
    print(f"[Self-Check] 合并样本数: {len(final_df)}")
    print(f"[Self-Check] 数据预览:\n{final_df.head(3)}")
    print("========================================================================")
 
 
if __name__ == "__main__":
    main()
