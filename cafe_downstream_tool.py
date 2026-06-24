#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import csv
from collections import defaultdict


# ---------- 工具函数 ----------

def count_nonempty_noncomment_lines(fname):
    """统计文件中非空且非注释行数。"""
    if not os.path.exists(fname):
        return 0
    n = 0
    with open(fname, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            n += 1
    return n


def read_family_set(fname):
    """读取家族ID文件（取第一列），返回set."""
    s = set()
    if not os.path.exists(fname):
        return s
    with open(fname, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fam = line.split("\t")[0].strip()
            if fam:
                s.add(fam)
    return s


# ---------- 模式 1：summary ----------

def run_mode_summary(species_list, cafe_dir, output):
    """
    功能1：
      自动生成汇总每个物种的扩张/收缩/显著扩张/显著收缩数量的表格 (TSV)。
    依赖文件（由 cafe_multi_tool.py sig 模式产生）：
      <cafe_dir>/<sp>_expand.txt
      <cafe_dir>/<sp>_contract.txt
      <cafe_dir>/<sp>_SigExpand.txt
      <cafe_dir>/<sp>_SigContract.txt
    """
    with open(output, "w", encoding="utf-8") as out:
        out.write("Species\tExpanded\tContracted\tSigExpanded\tSigContracted\n")
        for sp in species_list:
            f_exp = os.path.join(cafe_dir, f"{sp}_expand.txt")
            f_con = os.path.join(cafe_dir, f"{sp}_contract.txt")
            f_sigexp = os.path.join(cafe_dir, f"{sp}_SigExpand.txt")
            f_sigcon = os.path.join(cafe_dir, f"{sp}_SigContract.txt")

            n_exp = count_nonempty_noncomment_lines(f_exp)
            n_con = count_nonempty_noncomment_lines(f_con)
            n_sigexp = count_nonempty_noncomment_lines(f_sigexp)
            n_sigcon = count_nonempty_noncomment_lines(f_sigcon)

            if n_exp == 0 and n_con == 0 and n_sigexp == 0 and n_sigcon == 0:
                sys.stderr.write(
                    f"Warning: no data found for species {sp} "
                    f"({f_exp}, {f_con}, {f_sigexp}, {f_sigcon}).\n"
                )

            out.write(f"{sp}\t{n_exp}\t{n_con}\t{n_sigexp}\t{n_sigcon}\n")


# ---------- 模式 2：matrix ----------

def run_mode_matrix(species_list, cafe_dir, prefix="cafe_matrix"):
    """
    功能2：
      构建 “家族 × 物种” 的 0/1 矩阵，
      分别输出：
        prefix_expand.tsv
        prefix_contract.tsv
        prefix_SigExpand.tsv
        prefix_SigContract.tsv

    每个矩阵格式：
      FamilyID \t sp1 \t sp2 \t ...
    """
    categories = ["expand", "contract", "SigExpand", "SigContract"]

    # 为每个物种、每个类别读取集合
    fam_sets = {cat: {} for cat in categories}
    for sp in species_list:
        for cat in categories:
            fname = os.path.join(cafe_dir, f"{sp}_{cat}.txt")
            s = read_family_set(fname)
            if not s:
                sys.stderr.write(
                    f"Warning: file {fname} not found or empty for species {sp}.\n"
                )
            fam_sets[cat][sp] = s

    # 对每个类别构建家族全集并输出矩阵
    for cat in categories:
        all_fams = set()
        for sp in species_list:
            all_fams |= fam_sets[cat].get(sp, set())
        all_fams = sorted(all_fams)
        if not all_fams:
            sys.stderr.write(f"Note: no families found for category {cat}.\n")

        outname = f"{prefix}_{cat}.tsv"
        with open(outname, "w", encoding="utf-8") as out:
            # header
            out.write("FamilyID\t" + "\t".join(species_list) + "\n")
            for fam in all_fams:
                row = [fam]
                for sp in species_list:
                    row.append("1" if fam in fam_sets[cat].get(sp, set()) else "0")
                out.write("\t".join(row) + "\n")


# ---------- 模式 3：shared_genes ----------

def run_mode_shared_genes(species_list, cafe_dir, orthofinder_file, focal_species=None):
    """
    功能3：
      根据提供的物种id，提取这些物种共享的扩张/收缩/显著扩张/显著收缩家族的基因列表。
      只输出“对应的基因ID”，一行一个，没有额外信息。

      例如：-s Aart Mmic
        → 计算 Aart 和 Mmic 共享的 expand / contract / SigExpand / SigContract 家族
        → 对每一类输出 Aart 的基因ID列表（去重），一行一个。

      focal_species:
        若为 None，则默认 species_list[0] 为焦点物种；
        若指定，则必须在 species_list 中。
    """

    if len(species_list) < 2:
        sys.stderr.write("Error: shared_genes mode requires at least 2 species.\n")
        return

    if focal_species is None:
        focal_species = species_list[0]
    if focal_species not in species_list:
        sys.stderr.write(
            f"Error: focal species {focal_species} must be in the species list.\n"
        )
        return

    categories = ["expand", "contract", "SigExpand", "SigContract"]
    fam_sets = {cat: {} for cat in categories}

    # 1. 读取每个物种的家族集合（来自 cafe_multi_tool 的结果目录）
    for sp in species_list:
        for cat in categories:
            fname = os.path.join(cafe_dir, f"{sp}_{cat}.txt")
            s = read_family_set(fname)
            if not s:
                sys.stderr.write(
                    f"Warning: file {fname} not found or empty for species {sp}.\n"
                )
            fam_sets[cat][sp] = s

    # 2. 对每个类别计算“所有物种的交集家族”
    shared_fams = {}
    for cat in categories:
        sets_for_cat = [fam_sets[cat].get(sp, set()) for sp in species_list]
        if not sets_for_cat:
            shared_fams[cat] = set()
            continue
        inter = sets_for_cat[0].copy()
        for s in sets_for_cat[1:]:
            inter &= s
        shared_fams[cat] = inter
        sys.stderr.write(
            f"shared_genes: category {cat}, shared families among "
            f"{','.join(species_list)} = {len(inter)}\n"
        )

    # 3. 读取 Orthogroups.tsv，提取焦点物种的基因表
    with open(orthofinder_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            sys.stderr.write("Error: Orthogroups.tsv is empty.\n")
            return
        if focal_species not in header:
            sys.stderr.write(
                f"Error: focal species {focal_species} not found in Orthogroups header.\n"
            )
            return
        sp_idx = header.index(focal_species)

        og2genes = {}
        for row in reader:
            if not row:
                continue
            og = row[0].strip()
            if not og:
                continue
            if sp_idx >= len(row):
                cell = ""
            else:
                cell = row[sp_idx].strip()
            og2genes[og] = cell  # e.g. "gene1, gene2, gene3" or ""

    # 4. 对每个类别，收集所有共享家族中焦点物种的基因ID → 去重 → 每行一个ID
    others = [sp for sp in species_list if sp != focal_species]
    others_tag = "_".join(others) if others else "none"

    for cat in categories:
        fams = sorted(shared_fams.get(cat, set()))
        gene_ids = set()
        for fam in fams:
            genes_str = og2genes.get(fam, "")
            if not genes_str:
                continue
            # Orthofinder 中使用 ", " 作为分隔符
            for g in genes_str.split(","):
                g = g.strip()
                if g:
                    gene_ids.add(g)

        outname = f"{focal_species}_shared_with_{others_tag}_{cat}_genes.txt"
        with open(outname, "w", encoding="utf-8") as out:
            # 只输出对应ID，一行一个
            for gid in sorted(gene_ids):
                out.write(gid + "\n")


# ---------- 模式 4：idmap（GFF ID 转换） ----------

def parse_gff_attributes(attr_str):
    """将 GFF 第9列解析为 dict。"""
    attrs = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def build_id_maps_from_gff(gff_file):
    """
    根据 GFF 构建：
      gene_id -> mRNA_ids
      mRNA_id -> gene_id
      mRNA_id -> protein_ids
      protein_id -> mRNA_ids

    这里的 protein_id 优先使用 attributes 里的 'protein_id'，
    若没有，则在 CDS/polypeptide/protein feature 中用 ID 作为 protein_id 兜底。
    """
    gene_to_mrna = defaultdict(set)
    mrna_to_gene = {}
    mrna_to_protein = defaultdict(set)
    protein_to_mrna = defaultdict(set)

    gene_ids = set()
    mrna_ids = set()

    # 第一遍扫描：识别 gene / mRNA 以及 mRNA->gene
    with open(gff_file, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            attrs = parse_gff_attributes(parts[8])

            if ftype.lower() == "gene":
                gid = attrs.get("ID")
                if gid:
                    gene_ids.add(gid)
            elif ftype.lower() in ("mrna", "transcript"):
                mid = attrs.get("ID")
                parent = attrs.get("Parent")
                if mid:
                    mrna_ids.add(mid)
                if mid and parent:
                    mrna_to_gene[mid] = parent
                    gene_to_mrna[parent].add(mid)

    # 第二遍扫描：识别 protein 相关信息（CDS/polypeptide/protein）
    with open(gff_file, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            attrs = parse_gff_attributes(parts[8])

            if ftype.lower() in ("cds", "polypeptide", "protein"):
                parent = attrs.get("Parent")
                if not parent:
                    continue
                protein_id = (
                    attrs.get("protein_id")
                    or attrs.get("proteinId")
                    or attrs.get("protein")
                    or attrs.get("ID")
                )
                if not protein_id:
                    continue
                mrna_to_protein[parent].add(protein_id)
                protein_to_mrna[protein_id].add(parent)

    return gene_to_mrna, mrna_to_gene, mrna_to_protein, protein_to_mrna


def run_mode_idmap(gff_file, map_type, input_file, output_file):
    """
    功能4：
      根据输入的gff文件，实现给定的文件内容中基因id、转录本id和蛋白id的相互转化。
      输出只包含映射得到的ID，一行一个，无额外内容（所有输入的映射结果取并集）。

    支持的 map_type:
      gene2mrna
      gene2pep
      mrna2gene
      mrna2pep
      pep2gene
      pep2mrna
    """
    valid_types = {
        "gene2mrna",
        "gene2pep",
        "mrna2gene",
        "mrna2pep",
        "pep2gene",
        "pep2mrna",
    }
    if map_type not in valid_types:
        sys.stderr.write(
            "Error: map_type must be one of: "
            + ", ".join(sorted(valid_types))
            + "\n"
        )
        return

    gene_to_mrna, mrna_to_gene, mrna_to_protein, protein_to_mrna = \
        build_id_maps_from_gff(gff_file)

    # 为 gene2pep / pep2gene 衍生映射
    gene_to_protein = defaultdict(set)
    for gid, mrnas in gene_to_mrna.items():
        for mid in mrnas:
            for pep in mrna_to_protein.get(mid, []):
                gene_to_protein[gid].add(pep)

    protein_to_gene = defaultdict(set)
    for pep, mrnas in protein_to_mrna.items():
        for mid in mrnas:
            gid = mrna_to_gene.get(mid)
            if gid:
                protein_to_gene[pep].add(gid)

    # 读取输入 ID 列表
    ids = []
    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.append(line)

    # 统一收集所有映射结果的并集
    all_mapped = set()
    for x in ids:
        mapped = set()
        if map_type == "gene2mrna":
            mapped |= gene_to_mrna.get(x, set())
        elif map_type == "gene2pep":
            mapped |= gene_to_protein.get(x, set())
        elif map_type == "mrna2gene":
            gid = mrna_to_gene.get(x)
            if gid:
                mapped.add(gid)
        elif map_type == "mrna2pep":
            mapped |= mrna_to_protein.get(x, set())
        elif map_type == "pep2gene":
            mapped |= protein_to_gene.get(x, set())
        elif map_type == "pep2mrna":
            mapped |= protein_to_mrna.get(x, set())

        all_mapped |= mapped

    # 输出：只包含对应的ID，一行一个
    with open(output_file, "w", encoding="utf-8") as out:
        for mid in sorted(all_mapped):
            out.write(mid + "\n")


# ---------- 参数解析 & 主函数 ----------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Downstream analysis tool for cafe_multi_tool.py results.\n"
            "Modes:\n"
            "  summary      - summarize per-species counts of expand/contract/SigExpand/SigContract\n"
            "  matrix       - build 0/1 family×species matrices for expand/contract/SigExpand/SigContract\n"
            "  shared_genes - extract genes in shared families (expand/contract/SigExpand/SigContract)\n"
            "  idmap        - convert between gene/mRNA/protein IDs based on a GFF file\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-m", "--mode", action="append",
        choices=["summary", "matrix", "shared_genes", "idmap"],
        help="Mode(s) to run. Can be given multiple times."
    )

    # cafe_multi_tool 结果目录
    parser.add_argument(
        "--cafe-dir", default=".",
        help="Directory containing cafe_multi_tool.py output files "
             "(e.g. Species_expand.txt etc.). Default: current directory."
    )

    # 共用参数
    parser.add_argument(
        "-s", "--species", nargs="+",
        help="Species list (for 'summary' & 'matrix' & 'shared_genes' modes)."
    )

    # summary 输出
    parser.add_argument(
        "--summary-output", default="cafe_summary_per_species.tsv",
        help="Output TSV for summary mode. Default: cafe_summary_per_species.tsv"
    )

    # matrix 输出前缀
    parser.add_argument(
        "--matrix-prefix", default="cafe_matrix",
        help="Prefix for matrix mode output files. Default: cafe_matrix"
    )

    # shared_genes 相关
    parser.add_argument(
        "-of", "--orthofinder",
        help="Orthogroups.tsv from OrthoFinder (for 'shared_genes' mode)."
    )
    parser.add_argument(
        "--focal-species",
        help="Focal species for 'shared_genes' mode. If omitted, use the first species in -s."
    )

    # idmap 相关
    parser.add_argument(
        "--gff",
        help="GFF file for 'idmap' mode."
    )
    parser.add_argument(
        "--map-type",
        help="Mapping type for 'idmap' mode. "
             "One of: gene2mrna, gene2pep, mrna2gene, mrna2pep, pep2gene, pep2mrna."
    )
    parser.add_argument(
        "--id-input",
        help="Input file containing one ID per line (for 'idmap' mode)."
    )
    parser.add_argument(
        "--id-output", default="idmap_output.txt",
        help="Output file for 'idmap' mode. Each line: one mapped ID. Default: idmap_output.txt"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    modes = args.mode or []
    if not modes:
        sys.stderr.write("Error: no mode specified. Use -m/--mode.\n")
        sys.exit(1)

    cafe_dir = args.cafe_dir

    # summary
    if "summary" in modes:
        if not args.species:
            sys.stderr.write("Error: 'summary' mode requires --species.\n")
        else:
            run_mode_summary(args.species, cafe_dir, args.summary_output)

    # matrix
    if "matrix" in modes:
        if not args.species:
            sys.stderr.write("Error: 'matrix' mode requires --species.\n")
        else:
            run_mode_matrix(args.species, cafe_dir, prefix=args.matrix_prefix)

    # shared_genes
    if "shared_genes" in modes:
        if not args.species:
            sys.stderr.write("Error: 'shared_genes' mode requires --species.\n")
        elif not args.orthofinder:
            sys.stderr.write("Error: 'shared_genes' mode requires --orthofinder.\n")
        else:
            run_mode_shared_genes(
                args.species,
                cafe_dir,
                args.orthofinder,
                focal_species=args.focal_species
            )

    # idmap
    if "idmap" in modes:
        if not args.gff or not args.map_type or not args.id_input:
            sys.stderr.write(
                "Error: 'idmap' mode requires --gff, --map-type, and --id-input.\n"
            )
        else:
            run_mode_idmap(
                args.gff,
                args.map_type,
                args.id_input,
                args.id_output
            )


if __name__ == "__main__":
    main()

