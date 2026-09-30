import os
import argparse
import pandas as pd
import numpy as np
 
def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter TPM matrix by transcript-to-gene mapping and population-level expression threshold."
    )
    parser.add_argument("-m", "--tpm_matrix", required=True, help="Path to input salmon_tpm_matrix.csv")
    parser.add_argument("-pop", "--pop_info", required=True, help="Path to sample population mapping file (e.g., sample_pop_Lat.csv)")
    parser.add_argument("-p", "--psba_file", default="psba_id.txt", help="Path to psbA gene ID list file")
    parser.add_argument("-n", "--nlr_file", default="Aart_nlr_id.txt", help="Path to NLR gene ID list file")
    parser.add_argument(
        "--psba-pos",
        default=None,
        help="Path to positively selected psbA gene ID list file (optional, 1 gene ID per line)."
    )
    parser.add_argument("-o", "--output_dir", default="./output", help="Output directory path")
    return parser.parse_args()
 
def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 1. 读取原始基因 ID 列表
    psba_raw = [line.strip() for line in open(args.psba_file) if line.strip()]
    nlr_raw = [line.strip() for line in open(args.nlr_file) if line.strip()]
 
    # 处理可选的“受正选择 psbA” ID 列表
    pos_psba_genes = set()
    if args.psba_pos and os.path.exists(args.psba_pos):
        pos_psba_genes = set(line.strip() for line in open(args.psba_pos) if line.strip())
 
    raw_to_group = {gid: "psbA" for gid in psba_raw}
    raw_to_group.update({gid: "NLR" for gid in nlr_raw})
 
    # 2. 构建转录本 ID -> 基因 ID 的映射字典
    tx_to_gene = {}
    
    # 2.1 基础转换规则
    for gid in psba_raw + nlr_raw:
        if "evm.TU." in gid:
            tx_id = gid.replace("evm.TU.", "evm.model.")
        elif "novel_gene_" in gid:
            tx_id = gid.replace("novel_gene_", "novel_model_")
        else:
            tx_id = gid
        tx_to_gene[tx_id] = gid
        
    # 2.2 4 个基因的精确硬编码补丁 (修正带特定 Hash 后缀的转录本 ID)
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
 
    # 3. 读取已生成的表达矩阵与种群映射表
    tpm_matrix = pd.read_csv(args.tpm_matrix, index_col=0)
    
    pop_df = pd.read_csv(args.pop_info, sep=r'\s+|,', engine='python')
    sample_to_pop = dict(zip(pop_df['Sample'], pop_df['Population']))
 
    # 4. 执行 ID 匹配与转换
    matched_txs = [tx for tx in tx_to_gene.keys() if tx in tpm_matrix.index]
    matched_raw_ids = set(tx_to_gene[tx] for tx in matched_txs)
    
    unmatched_psba = set(psba_raw) - matched_raw_ids
    unmatched_nlr = set(nlr_raw) - matched_raw_ids
 
    # ================== Self-Check 检查 ==================
    print("========================================================================")
    print(f"[Self-Check] 构建映射规则数: {len(tx_to_gene)}")
    print(f"[Self-Check] 表达矩阵中精准匹配到的转录本数: {len(matched_txs)}")
    if unmatched_psba:
        print(f"[WARNING] 缺失 psbA 基因 ({len(unmatched_psba)} 个): {list(unmatched_psba)}")
    if unmatched_nlr:
        print(f"[WARNING] 缺失 NLR 基因 ({len(unmatched_nlr)} 个): {list(unmatched_nlr)}")
    if not unmatched_psba and not unmatched_nlr:
        print("[SUCCESS] 目标基因 100% 精确匹配！")
    print("========================================================================")
    # ======================================================
 
    assert len(matched_txs) > 0, "ERROR: 未在表达矩阵中匹配到任何目标基因，请检查 ID 列表或矩阵文件"
 
    # 提取目标基因子矩阵，将行名替换回真实的 Gene_ID
    sub_matrix = tpm_matrix.loc[matched_txs].copy()
    sub_matrix.index = [tx_to_gene[tx] for tx in sub_matrix.index]
    sub_matrix.index.name = "Gene_ID"
 
    # 5. 核心过滤逻辑：仅保留在至少 1 个种群中，存在 >= 2 个样本个体 TPM >= 1.0 的基因
    valid_samples = [col for col in sub_matrix.columns if col in sample_to_pop]
    sub_matrix = sub_matrix[valid_samples]
    pop_series = pd.Series({s: sample_to_pop[s] for s in valid_samples})
 
    bool_matrix = (sub_matrix >= 1.0)
    # 消除 Pandas 未来版本弃用 warning (采用矩阵转置置换)
    pop_counts = bool_matrix.T.groupby(pop_series).sum().T
 
    passed_mask = (pop_counts >= 2).any(axis=1)
    filtered_matrix = sub_matrix[passed_mask].copy()
 
    # 6. 生成 Meta 信息与输出
    gene_info = []
    for gene_id in filtered_matrix.index:
        gene_info.append({"Gene_ID": gene_id, "Group": raw_to_group[gene_id]})
    df_info = pd.DataFrame(gene_info)
 
    filtered_matrix.to_csv(os.path.join(args.output_dir, "filtered_tpm_matrix.csv"))
    df_info.to_csv(os.path.join(args.output_dir, "gene_info.csv"), index=False)
 
    long_df = filtered_matrix.reset_index().melt(
        id_vars="Gene_ID",
        var_name="Sample",
        value_name="TPM"
    )
    long_df = long_df.merge(df_info, on="Gene_ID")
    long_df["log2_TPM"] = np.log2(long_df["TPM"] + 1)
    long_df.to_csv(os.path.join(args.output_dir, "long_format.csv"), index=False)
 
    # 7. 打印汇总信息与 psbA 明细统计
    passed_psba = set(df_info[df_info['Group'] == 'psbA']['Gene_ID'])
    passed_nlr = set(df_info[df_info['Group'] == 'NLR']['Gene_ID'])
    
    print(f"[SUCCESS] 过滤后保留基因总数: {len(filtered_matrix)}")
    print(f"其中 psbA: {len(passed_psba)} 个, NLR: {len(passed_nlr)} 个")
 
    if args.psba_pos:
        passed_pos_psba = passed_psba.intersection(pos_psba_genes)
        passed_nonpos_psba = passed_psba - pos_psba_genes
        print(f"[INFO] psbA 分组明细:")
        print(f"       - 受正选择 psbA: 保留 {len(passed_pos_psba)} / {len(pos_psba_genes)} 个")
        print(f"       - 未受正选择 psbA: 保留 {len(passed_nonpos_psba)} 个")
 
if __name__ == "__main__":
    main()
