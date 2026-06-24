#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import sys
import os
import re
from collections import defaultdict


# ---------- 通用小工具 ----------

def safe_int(x):
    x = x.strip()
    if x == "" or x == ".":
        return 0
    try:
        return int(x)
    except ValueError:
        try:
            return int(float(x))
        except ValueError:
            raise ValueError(f"Cannot convert value to int: {x!r}")


# ---------- 解析 Gamma_change.tab ----------

def parse_gamma_change(gamma_change_file):
    """
    解析 Gamma_change.tab：
    返回：
      species_list: [sp1, sp2, ...]  （按列顺序）
      family_ids: [OG0000001, OG0000002, ...]
      values: {family: {sp: int}}
      expanded_fams: {sp: set(family)}
      contracted_fams: {sp: set(family)}
    """
    species_list = []
    values = {}
    expanded_fams = defaultdict(set)
    contracted_fams = defaultdict(set)
    family_ids = []

    with open(gamma_change_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            sys.stderr.write("Gamma_change.tab is empty.\n")
            sys.exit(1)

        # 从第二列开始识别物种列
        species_cols = []  # [(species, idx), ...]
        for idx, colname in enumerate(header[1:], start=1):
            colname = colname.strip()
            if not colname:
                continue
            if "<" in colname:
                sp = colname.split("<", 1)[0].strip()
            else:
                sp = colname.strip()
            if sp == "":
                # 跳过 "<17>" 这类没有物种名的列
                continue
            species_cols.append((sp, idx))
            species_list.append(sp)

        if not species_cols:
            sys.stderr.write("No species columns detected in Gamma_change header.\n")
            sys.exit(1)

        for row in reader:
            if not row:
                continue
            fam = row[0].strip()
            if not fam:
                continue
            family_ids.append(fam)
            values[fam] = {}

            for sp, idx in species_cols:
                if idx >= len(row):
                    v = 0
                else:
                    v = safe_int(row[idx])
                values[fam][sp] = v
                if v > 0:
                    expanded_fams[sp].add(fam)
                elif v < 0:
                    contracted_fams[sp].add(fam)
                # v == 0 不记录

    return species_list, family_ids, values, expanded_fams, contracted_fams


# ---------- 解析 Gamma_asr.tre ----------

def parse_gamma_asr(gamma_asr_file, species_list):
    """
    只看叶节点（物种），判定是否显著（带 *）。
    叶节点形如：
        Ccan<1>_5:...
        Alap<2>*_21:...

    对每个家族、每个物种，检测是否匹配：
        species<digits>(*?)_
    若 group(1) == "*" 则该家族在该物种中显著变化。

    返回：
      sig_fams: {sp: set(family_id_with_star)}
    """
    sig_fams = {sp: set() for sp in species_list}

    # 预编译每个物种的正则
    species_patterns = {
        sp: re.compile(rf'{re.escape(sp)}<\d+>(\*?)_')
        for sp in species_list
    }

    with open(gamma_asr_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or not line.startswith("TREE"):
                continue
            # 解析 TREE 名和树串
            # 例: TREE OG0000063 = ( ... );
            m = re.match(r'TREE\s+(\S+)\s*=\s*(.*);', line, flags=re.IGNORECASE)
            if not m:
                continue
            fam = m.group(1)
            tree_str = m.group(2)

            for sp, pattern in species_patterns.items():
                m2 = pattern.search(tree_str)
                if m2 and m2.group(1) == "*":
                    sig_fams[sp].add(fam)

    return sig_fams


# ---------- 模式 1：对所有物种输出 expand/contract + SigExpand/SigContract ----------

def run_mode_sig(species_list, expanded_fams, contracted_fams, sig_fams):
    """
    对所有物种：
      1) 输出扩张家族: 物种_expand.txt
      2) 输出收缩家族: 物种_contract.txt
      3) 输出显著扩张家族: 物种_SigExpand.txt
      4) 输出显著收缩家族: 物种_SigContract.txt
    """
    for sp in species_list:
        exp_set = expanded_fams.get(sp, set())
        con_set = contracted_fams.get(sp, set())
        sig_set = sig_fams.get(sp, set())

        sig_exp_set = exp_set & sig_set
        sig_con_set = con_set & sig_set

        # 普通扩张
        with open(f"{sp}_expand.txt", "w", encoding="utf-8") as out:
            for fam in sorted(exp_set):
                out.write(f"{fam}\n")

        # 普通收缩
        with open(f"{sp}_contract.txt", "w", encoding="utf-8") as out:
            for fam in sorted(con_set):
                out.write(f"{fam}\n")

        # 显著扩张
        with open(f"{sp}_SigExpand.txt", "w", encoding="utf-8") as out:
            for fam in sorted(sig_exp_set):
                out.write(f"{fam}\n")

        # 显著收缩
        with open(f"{sp}_SigContract.txt", "w", encoding="utf-8") as out:
            for fam in sorted(sig_con_set):
                out.write(f"{fam}\n")


# ---------- 模式 2：指定物种集合的共享 expand/contract + SigExpand/SigContract ----------

def run_mode_shared(requested_species, species_list, family_ids,
                    expanded_fams, contracted_fams, sig_fams):
    """
    requested_species: 用户通过 --species 指定的物种列表
    species_list: Gamma_change 中真实存在的物种列表
    family_ids: 所有家族 ID 列表
    expanded_fams, contracted_fams: {sp: set(fams)}
    sig_fams: {sp: set(fams_with_star)} (ASR)

    输出：
      Nsp_expand.txt / Nsp_contract.txt
      Nsp_SigExpand.txt / Nsp_SigContract.txt
    格式：
      对长度为 N 的共享家族：只输出 FamilyID
      对长度为 2..N-1 的共享家族：FamilyID \t sp1+sp2+...
    """
    # 确认哪些 species 有效
    valid_species = [sp for sp in requested_species if sp in species_list]
    invalid_species = [sp for sp in requested_species if sp not in species_list]

    if invalid_species:
        sys.stderr.write(
            "Warning: the following species not found in Gamma_change and will be ignored: "
            + ", ".join(invalid_species) + "\n"
        )

    if len(valid_species) < 2:
        sys.stderr.write("Need at least 2 valid species for shared analysis.\n")
        return

    N = len(valid_species)

    # 显著扩张/收缩 per species
    sig_expanded_fams = {sp: expanded_fams.get(sp, set()) & sig_fams.get(sp, set())
                         for sp in valid_species}
    sig_contracted_fams = {sp: contracted_fams.get(sp, set()) & sig_fams.get(sp, set())
                           for sp in valid_species}

    shared_expand = defaultdict(list)       # k -> [(fam, [sp...])]
    shared_contract = defaultdict(list)
    shared_sig_expand = defaultdict(list)
    shared_sig_contract = defaultdict(list)

    for fam in family_ids:
        # 普通扩张
        ex_sp = [sp for sp in valid_species if fam in expanded_fams.get(sp, set())]
        if len(ex_sp) >= 2:
            shared_expand[len(ex_sp)].append((fam, ex_sp))

        # 普通收缩
        con_sp = [sp for sp in valid_species if fam in contracted_fams.get(sp, set())]
        if len(con_sp) >= 2:
            shared_contract[len(con_sp)].append((fam, con_sp))

        # 显著扩张
        sig_ex_sp = [sp for sp in valid_species if fam in sig_expanded_fams.get(sp, set())]
        if len(sig_ex_sp) >= 2:
            shared_sig_expand[len(sig_ex_sp)].append((fam, sig_ex_sp))

        # 显著收缩
        sig_con_sp = [sp for sp in valid_species if fam in sig_contracted_fams.get(sp, set())]
        if len(sig_con_sp) >= 2:
            shared_sig_contract[len(sig_con_sp)].append((fam, sig_con_sp))

    # 输出文件名
    expand_shared_file = f"{N}sp_expand.txt"
    contract_shared_file = f"{N}sp_contract.txt"
    sig_expand_shared_file = f"{N}sp_SigExpand.txt"
    sig_contract_shared_file = f"{N}sp_SigContract.txt"

    # 普通扩张
    with open(expand_shared_file, "w", encoding="utf-8") as out:
        for k in range(N, 1, -1):
            fam_list = shared_expand.get(k, [])
            if not fam_list:
                continue
            for fam, sp_list in sorted(fam_list, key=lambda x: x[0]):
                if k == N:
                    out.write(f"{fam}\n")
                else:
                    sp_str = "+".join(sp_list)
                    out.write(f"{fam}\t{sp_str}\n")

    # 普通收缩
    with open(contract_shared_file, "w", encoding="utf-8") as out:
        for k in range(N, 1, -1):
            fam_list = shared_contract.get(k, [])
            if not fam_list:
                continue
            for fam, sp_list in sorted(fam_list, key=lambda x: x[0]):
                if k == N:
                    out.write(f"{fam}\n")
                else:
                    sp_str = "+".join(sp_list)
                    out.write(f"{fam}\t{sp_str}\n")

    # 显著扩张
    with open(sig_expand_shared_file, "w", encoding="utf-8") as out:
        for k in range(N, 1, -1):
            fam_list = shared_sig_expand.get(k, [])
            if not fam_list:
                continue
            for fam, sp_list in sorted(fam_list, key=lambda x: x[0]):
                if k == N:
                    out.write(f"{fam}\n")
                else:
                    sp_str = "+".join(sp_list)
                    out.write(f"{fam}\t{sp_str}\n")

    # 显著收缩
    with open(sig_contract_shared_file, "w", encoding="utf-8") as out:
        for k in range(N, 1, -1):
            fam_list = shared_sig_contract.get(k, [])
            if not fam_list:
                continue
            for fam, sp_list in sorted(fam_list, key=lambda x: x[0]):
                if k == N:
                    out.write(f"{fam}\n")
                else:
                    sp_str = "+".join(sp_list)
                    out.write(f"{fam}\t{sp_str}\n")


# ---------- 模式 3：从 Orthogroups.tsv + 家族ID文件 提取基因 ----------

def autodetect_family_files_for_species(sp):
    """
    自动检测某个物种对应的四类家族文件：
      sp_expand.txt
      sp_contract.txt
      sp_SigExpand.txt
      sp_SigContract.txt
    返回存在的文件列表（按上述顺序）。
    """
    candidates = [
        f"{sp}_expand.txt",
        f"{sp}_contract.txt",
        f"{sp}_SigExpand.txt",
        f"{sp}_SigContract.txt",
    ]
    return [fn for fn in candidates if os.path.exists(fn)]


def run_mode_extract(orthofinder_file, species_list, family_files, num_extract_species):
    """
    species_list: 通过 -s 指定的物种列表
    family_files: 用户手动指定的家族ID文件列表（可为 None）
    num_extract_species: 通过 -n 指定的前 N 个物种

    Orthogroups.tsv:
      第一列是 Orthogroup ID
      其他列为物种名，对应的单元格是 "gene1, gene2, gene3"（逗号+空格）

    输出文件名：物种名_家族id文件名(去掉.txt).txt
      例如 species=Ccan, family_file=Ccan_expand.txt -> Ccan_Ccan_expand.txt

    输出格式（每行）：
      FamilyID \t gene1, gene2, gene3
    """
    if not species_list:
        sys.stderr.write("Error: 'extract' mode requires --species.\n")
        return

    # 实际要处理的物种集合：-s 的前 N 个，N 由 -n 决定，默认 1
    num = max(1, num_extract_species)
    num = min(num, len(species_list))
    extract_species = species_list[:num]

    # 先解析 Orthogroups.tsv
    with open(orthofinder_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            sys.stderr.write("Orthogroups.tsv is empty.\n")
            return

        # 对每个物种，找到其列索引
        sp2idx = {}
        for sp in extract_species:
            if sp not in header:
                sys.stderr.write(
                    f"Warning: Species {sp} not found in Orthogroups header. It will be skipped in extract mode.\n"
                )
            else:
                sp2idx[sp] = header.index(sp)

        if not sp2idx:
            sys.stderr.write("Error: None of the requested species were found in Orthogroups.tsv header.\n")
            return

        # 将 Orthogroups.tsv 全部读入内存（对大文件可能稍重，但通常可接受）
        og_rows = []
        for row in reader:
            if not row:
                continue
            og = row[0].strip()
            if not og:
                continue
            og_rows.append(row)

    # 为了避免多次遍历 og_rows，这里构建一个字典：og -> {sp: cell_str}
    og2genes_per_sp = {sp: {} for sp in sp2idx.keys()}

    for row in og_rows:
        og = row[0].strip()
        for sp, idx in sp2idx.items():
            if idx >= len(row):
                cell = ""
            else:
                cell = row[idx].strip()
            og2genes_per_sp[sp][og] = cell  # cell 如 "g1, g2, g3" 或 ""

    # 逐个物种、逐个家族文件处理
    for sp in extract_species:
        if sp not in sp2idx:
            continue  # 上面已经警告过

        # 确定该物种要用的家族文件列表
        if family_files:
            fam_files_for_sp = family_files
        else:
            fam_files_for_sp = autodetect_family_files_for_species(sp)
            if not fam_files_for_sp:
                sys.stderr.write(
                    f"Warning: No family files provided and no auto-detected files found for species {sp}.\n"
                )
                continue

        for fam_file in fam_files_for_sp:
            # 读取家族ID列表（取第一列）
            fam_ids = []
            seen = set()
            with open(fam_file, "r", encoding="utf-8") as fin:
                for line in fin:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    cols = line.split("\t")
                    fam = cols[0].strip()
                    if fam and fam not in seen:
                        seen.add(fam)
                        fam_ids.append(fam)

            base = os.path.basename(fam_file)
            if base.lower().endswith(".txt"):
                base = base[:-4]
            out_name = f"{sp}_{base}.txt"

            with open(out_name, "w", encoding="utf-8") as out:
                for fam in fam_ids:
                    genes_str = og2genes_per_sp[sp].get(fam, "")
                    # 保持 ", " 这种分隔符原样输出
                    out.write(f"{fam}\t{genes_str}\n")


# ---------- 主程序 & 参数解析 ----------
LONG_HELP = """
============================================================
         Extended Help for CAFE / ASR / OrthoFinder Tool
============================================================

This script provides an integrated workflow for analyzing
gene family expansion/contraction (CAFE Gamma_change),
significant ancestral shifts (CAFE Gamma_asr),
and extracting gene IDs from OrthoFinder results.

You may run modes independently using --mode, or run all 
modes automatically using --all (in the order: 
sig → shared → extract).

------------------------------------------------------------
1) Mode: sig  (per-species expansion / contraction summary)
------------------------------------------------------------

For every species detected in Gamma_change.tab, the script 
outputs:

  Species_expand.txt
      Families with positive change value (expansion)

  Species_contract.txt
      Families with negative change value (contraction)

  Species_SigExpand.txt
      Families that are BOTH:
          - expanded in Gamma_change.tab
          - marked with "*" as significant in Gamma_asr.tre

  Species_SigContract.txt
      Families that are BOTH:
          - contracted in Gamma_change.tab
          - marked with "*" in Gamma_asr.tre

Requirements:
  --gamma-change
  --gamma-asr

------------------------------------------------------------
2) Mode: shared  (shared expansions among selected species)
------------------------------------------------------------

Given multiple species via --species sp1 sp2 sp3 ...

Outputs:

  Nsp_expand.txt
  Nsp_contract.txt
  Nsp_SigExpand.txt
  Nsp_SigContract.txt

Where N is the number of valid species.

Output rule:
  - If a family is shared by ALL species: 
        print only FamilyID
  - If shared by k < N species:
        FamilyID    sp1+sp2+sp3

Requirements:
  --gamma-change
  --gamma-asr
  --species sp1 sp2 ...

------------------------------------------------------------
3) Mode: extract  (extract gene IDs from Orthogroups.tsv)
------------------------------------------------------------

This mode extracts gene IDs belonging to specific families.

Inputs:
  --orthofinder  Path to Orthogroups.tsv
  --species      One or more species
  --family-files Optional family ID lists

If no --family-files are provided, the script automatically 
uses the 4 files created by sig mode:

  Species_expand.txt
  Species_contract.txt
  Species_SigExpand.txt
  Species_SigContract.txt

Species selection:
  - Only the first species in -s is used by default
  - Use -n N to process the first N species instead

For species Sp and family file Sp_expand.txt:
  Output → Sp_Sp_expand.txt

Each line:
  FamilyID    gene1, gene2, gene3

Requirements:
  --orthofinder
  --species

------------------------------------------------------------
4) --all  (run everything automatically)
------------------------------------------------------------

Equivalent to running:

    --mode sig 
    --mode shared 
    --mode extract

All mode dependencies are handled automatically.

------------------------------------------------------------
Examples
------------------------------------------------------------

Run sig mode:
  python script.py -gc Gamma_change.tab -ga Gamma_asr.tre -m sig

Run shared mode:
  python script.py -gc GC.tab -ga ASR.tre -s Ccan Alap Csol -m shared

Run extract mode (first 2 species):
  python script.py -of Orthogroups.tsv -s Aart Cint Ecan -n 2 -m extract

Run all modes:
  python script.py -gc GC.tab -ga ASR.tre -of OG.tsv -s Ccan Alap Csol -all

============================================================
"""

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Multi-mode CAFE/ASR/Orthofinder helper.\n"
            "Use --help-long for detailed documentation."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-gc", "--gamma-change",
                        help="Gamma_change.tab from CAFE")
    parser.add_argument("-ga", "--gamma-asr",
                        help="Gamma_asr.tre from CAFE")
    parser.add_argument("-of", "--orthofinder",
                        help="Orthogroups.tsv from OrthoFinder")

    parser.add_argument("-s", "--species", nargs="+",
                        help="Species list used in shared/extract modes.")
    parser.add_argument("-ff", "--family-files", nargs="+",
                        help="Family ID files for extract mode.")
    parser.add_argument("-n", "--num-extract-species", type=int, default=1,
                        help="Number of species (from -s list) to process in extract mode. Default = 1.")

    parser.add_argument("-m", "--mode", action="append",
                        choices=["sig", "shared", "extract"],
                        help="Run specific mode(s): sig, shared, extract")
    parser.add_argument("-all", action="store_true",
                        help="Run all modes in order: sig → shared → extract")

    parser.add_argument("--help-long", action="store_true",
                        help="Show extended documentation")

    return parser.parse_args()




def main():


    args = parse_args()

    if args.help_long:
        print(LONG_HELP)
        sys.exit(0)
    # 决定要运行的模式列表
    if args.all:
        modes = ["sig", "shared", "extract"]
    else:
        modes = args.mode or []

    if not modes:
        sys.stderr.write("Error: No mode specified. Use --mode or --all.\n")
        sys.exit(1)

    # 需要 Gamma_change / Gamma_asr 的模式？
    need_change = any(m in ("sig", "shared") for m in modes)
    need_asr = any(m in ("sig", "shared") for m in modes)

    # 存储解析结果（避免重复解析）
    species_list = None
    family_ids = None
    values = None
    expanded_fams = None
    contracted_fams = None
    sig_fams = None

    # 解析 Gamma_change.tab
    if need_change:
        if not args.gamma_change:
            sys.stderr.write("Error: --gamma-change is required for 'sig' and 'shared' modes.\n")
            sys.exit(1)
        (species_list, family_ids, values,
         expanded_fams, contracted_fams) = parse_gamma_change(args.gamma_change)

    # 解析 Gamma_asr.tre
    if need_asr:
        if not args.gamma_asr:
            sys.stderr.write("Error: --gamma-asr is required for 'sig' and 'shared' modes.\n")
            sys.exit(1)
        sig_fams = parse_gamma_asr(args.gamma_asr, species_list)

    # 依次按 sig -> shared -> extract 顺序运行
    # 但只对用户实际选择的模式生效
    if "sig" in modes:
        run_mode_sig(species_list, expanded_fams, contracted_fams, sig_fams)

    if "shared" in modes:
        if not args.species:
            sys.stderr.write("Error: --species is required for 'shared' mode.\n")
        else:
            run_mode_shared(args.species, species_list, family_ids,
                            expanded_fams, contracted_fams, sig_fams)

    if "extract" in modes:
        if not args.orthofinder:
            sys.stderr.write("Error: --orthofinder is required for 'extract' mode.\n")
        elif not args.species:
            sys.stderr.write("Error: 'extract' mode requires --species.\n")
        else:
            run_mode_extract(args.orthofinder,
                             args.species,
                             args.family_files,
                             args.num_extract_species)


if __name__ == "__main__":
    main()

