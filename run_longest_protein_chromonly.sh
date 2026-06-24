#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------
# Usage:
#   bash run_longest_protein_chromonly.sh genome.fa annotation.gff3 out_prefix [wrap_width] [min_chr_len] [cum_cov]
#
# Example:
#   bash run_longest_protein_chromonly.sh Alap_genome.fa Alap_genome.gff Alap 80 1000000 0.9
#
# Parameters:
#   wrap_width  : FASTA折行宽度(默认80)
#   min_chr_len : 认为是染色体的最短长度(默认1,000,000 bp)
#   cum_cov     : 最长contig累计覆盖比例阈值(默认0.9)
#
# Outputs:
#   ${PREFIX}.w80.fa
#   ${PREFIX}.chrom.list.txt                # 自动识别出的染色体/伪染色体ID
#   ${PREFIX}.chromonly.gff3                # 只保留这些ID上的注释
#   ${PREFIX}.chromonly.longest.uniqCDS.gff3
#   ${PREFIX}.chromonly.cds.fa
#   ${PREFIX}.chromonly.protein.longest.fa
# -----------------------------------------

FASTA="$1"
GFF="$2"
PREFIX="$3"
WRAP_WIDTH="${4:-80}"
MIN_CHR_LEN="${5:-1000000}"
CUM_COV="${6:-0.9}"

# check deps
for cmd in seqkit samtools agat_sp_keep_longest_isoform.pl agat_sp_extract_sequences.pl python; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "[ERROR] missing command: $cmd"; exit 1; }
done

echo "== Step 1: Wrap FASTA to ${WRAP_WIDTH} bp/line =="
WRAPPED_FASTA="${PREFIX}.w${WRAP_WIDTH}.fa"
if [[ ! -s "$WRAPPED_FASTA" ]]; then
  seqkit seq -w "$WRAP_WIDTH" "$FASTA" > "$WRAPPED_FASTA"
else
  echo "  [skip] wrapped FASTA exists: $WRAPPED_FASTA"
fi

echo "== Step 2: Rebuild FASTA index (.fai) =="
rm -f "${WRAPPED_FASTA}.fai"
samtools faidx "$WRAPPED_FASTA"

echo "== Step 3: Infer chromosome/pseudochromosome contigs from FASTA =="
CHR_LIST="${PREFIX}.chrom.list.txt"

python - "$WRAPPED_FASTA.fai" "$CHR_LIST" "$MIN_CHR_LEN" "$CUM_COV" <<'PY'
import sys, re

fai, out, min_len, cum_cov = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])

contigs = []
total = 0
with open(fai) as f:
    for ln in f:
        c, l = ln.split()[:2]
        l = int(l)
        contigs.append((c, l))
        total += l

contigs.sort(key=lambda x: x[1], reverse=True)

# 1) 长度累计到 cum_cov
keep = set()
cum = 0
max_len = contigs[0][1] if contigs else 0

for c, l in contigs:
    if l < min_len:
        continue
    cum += l
    keep.add(c)
    if cum / total >= cum_cov:
        break

# 2) 命名兜底（常见染色体/伪染色体命名）
pat = re.compile(r"(chr|chromosome|linkage|lg\d+|hic_asm|pseudo)", re.I)
for c, l in contigs:
    if l >= min_len and pat.search(c):
        keep.add(c)

with open(out, "w") as fw:
    for c, l in contigs:
        if c in keep:
            fw.write(c + "\n")

sys.stderr.write(f"[INFO] Total assembly size={total}\n")
sys.stderr.write(f"[INFO] Kept contigs={len(keep)} (min_len={min_len}, cum_cov={cum_cov})\n")
PY

echo "  Chromosome-like contigs saved to: $CHR_LIST"
echo "  Preview:"
head "$CHR_LIST"

echo "== Step 4: Filter GFF to chromosome-like contigs only =="
CHROMONLY_GFF="${PREFIX}.chromonly.gff3"

python - "$GFF" "$CHR_LIST" "$CHROMONLY_GFF" <<'PY'
import sys

gff, chr_list, out = sys.argv[1], sys.argv[2], sys.argv[3]
keep = set(x.strip() for x in open(chr_list) if x.strip())

with open(gff) as fin, open(out, "w") as fout:
    for ln in fin:
        if ln.startswith("#") or ln.strip()=="":
            fout.write(ln); continue
        cols = ln.split("\t")
        if cols[0] in keep:
            fout.write(ln)
PY

echo "== Step 5: Keep longest isoform per gene (AGAT) on chrom-only GFF =="
LONGEST_GFF="${PREFIX}.chromonly.longest.gff3"
agat_sp_keep_longest_isoform.pl -gff "$CHROMONLY_GFF" -o "$LONGEST_GFF"

echo "== Step 6: Make CDS IDs unique (Python) =="
UNIQ_GFF="${PREFIX}.chromonly.longest.uniqCDS.gff3"

python - "$LONGEST_GFF" "$UNIQ_GFF" <<'PY'
import sys, re
from collections import defaultdict

gff_in, gff_out = sys.argv[1], sys.argv[2]
counter = defaultdict(int)

def get_parent(attr):
    m = re.search(r"Parent=([^;]+)", attr)
    return m.group(1) if m else "NA"

def get_id(attr):
    m = re.search(r"ID=([^;]+)", attr)
    return m.group(1) if m else None

with open(gff_in) as fin, open(gff_out, "w") as fout:
    for line in fin:
        if line.startswith("#") or line.strip()=="":
            fout.write(line); continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            fout.write(line); continue

        if cols[2] == "CDS":
            attr = cols[8]
            parent = get_parent(attr)
            counter[parent] += 1
            idx = counter[parent]

            base_id = get_id(attr) or f"cds.{parent}"
            new_id = f"{base_id}.{idx}"

            if "ID=" in attr:
                attr = re.sub(r"ID=[^;]+", f"ID={new_id}", attr)
            else:
                attr = f"ID={new_id};" + attr

            cols[8] = attr
            fout.write("\t".join(cols) + "\n")
        else:
            fout.write(line)
PY

echo "== Step 7: Extract spliced CDS (AGAT) from chrom-only longest GFF =="
CDS_FA="${PREFIX}.chromonly.cds.fa"
agat_sp_extract_sequences.pl \
  --gff "$UNIQ_GFF" \
  --fasta "$WRAPPED_FASTA" \
  --type CDS \
  -o "$CDS_FA"

echo "== Step 8: Translate CDS -> protein (seqkit) =="
PROT_FA="${PREFIX}.chromonly.protein.longest.fa"
seqkit translate "$CDS_FA" > "$PROT_FA"

echo "== DONE =="
echo "Wrapped FASTA : $WRAPPED_FASTA"
echo "Chrom contigs : $CHR_LIST"
echo "Chrom-only GFF: $CHROMONLY_GFF"
echo "Longest GFF   : $LONGEST_GFF"
echo "Unique CDS GFF: $UNIQ_GFF"
echo "CDS FASTA     : $CDS_FA"
echo "Protein FASTA : $PROT_FA"

