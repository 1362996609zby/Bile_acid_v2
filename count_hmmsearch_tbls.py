# conda activate hmm_env
# python count_hmmsearch_tbls_filtered.py

import os
import pandas as pd

# ========= 用户参数 =========
root_dir = os.getcwd()
EVALUE_THRESHOLD = 1e-3   # <<< 关键修正：与 prior-check 一致
# ============================

results = {}

for dirpath, dirnames, filenames in os.walk(root_dir):
    for fn in filenames:
        if not fn.endswith(".tbl"):
            continue

        tbl_path = os.path.join(dirpath, fn)

        # 第一层目录作为 gene_dir
        rel = os.path.relpath(tbl_path, root_dir)
        parts = rel.split(os.sep)
        if len(parts) < 2:
            continue
        gene_dir = parts[0]

        bin_name = os.path.splitext(fn)[0]

        if bin_name not in results:
            results[bin_name] = {}

        hit_count = 0

        with open(tbl_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split()
                if len(fields) < 6:
                    continue  # 非标准 tbl 行跳过

                try:
                    full_evalue = float(fields[4])  # HMMER tblout 第5列 = full-seq E-value
                except:
                    continue

                if full_evalue <= EVALUE_THRESHOLD:
                    hit_count += 1

        results[bin_name][gene_dir] = results[bin_name].get(gene_dir, 0) + hit_count

df = pd.DataFrame.from_dict(results, orient="index").fillna(0).astype(int)

outname = f"hmmsearch_counts_E{EVALUE_THRESHOLD:.0e}.xlsx"
df.to_excel(outname)

print(f"✅ 统计完成（E-value ≤ {EVALUE_THRESHOLD}），结果保存在 {outname}")
