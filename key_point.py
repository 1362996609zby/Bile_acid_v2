# conda activate cdhit_env
# python key_point.py
# cat BaiN_ncbi.fasta BaiN_uniprot.fasta > BaiN_homo.fasta
import json
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# **1️⃣ 解析 *.json，提取关键位点**
def parse_json(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)

    key_sites = []
    for feature in data["features"]:
        if feature["type"] in ["Active site", "Binding site"]:
            start = feature["location"]["start"]["value"]
            end = feature["location"]["end"]["value"]
            for pos in range(start, end + 1):
                key_sites.append(pos)

    return sorted(set(key_sites))

# **2️⃣ 确保核心序列 *_core 始终存在**
def ensure_core_sequence(core_fasta, output_fasta):
    core_seq_record = list(SeqIO.parse(core_fasta, "fasta"))[0]  # 读取核心序列

    # 读取当前文件的所有序列
    records = list(SeqIO.parse(output_fasta, "fasta"))
    record_ids = [record.id for record in records]

    # 如果核心序列不在结果中，则添加
    if core_seq_record.id not in record_ids:
        print(f"⚠ 强制保留核心序列 {core_seq_record.id}")
        records.append(core_seq_record)

    # 重新写入文件，确保核心序列保留
    with open(output_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")

# **3️⃣ 运行 BLAST**
def run_blast(query_fasta, db_fasta, output_file, output_dir):
    db_name = os.path.join(output_dir, "blast_db")
    
    subprocess.run(["makeblastdb", "-in", db_fasta, "-dbtype", "prot", "-out", db_name])
    
    subprocess.run([
        "blastp", "-query", query_fasta, "-db", db_name, 
        "-out", output_file, "-evalue", "1e-5", "-max_target_seqs", "500", "-outfmt", "6"
    ])

# **4️⃣ 过滤 BLAST 结果，并检查关键位点，同时强制保留核心序列**
def filter_blast_results(blast_output, db_fasta, key_sites, filtered_fasta, core_fasta):
    valid_ids = set()
    
    with open(blast_output, "r") as f:
        for line in f:
            parts = line.split("\t")
            seq_id = parts[1]
            valid_ids.add(seq_id)

    with open(filtered_fasta, "w") as f_out:
        for record in SeqIO.parse(db_fasta, "fasta"):
            if record.id in valid_ids:
                seq = str(record.seq).upper()
                if all(seq[pos - 1].isupper() for pos in key_sites if pos - 1 < len(seq)):
                    SeqIO.write(record, f_out, "fasta")
                else:
                    print(f"⚠ 过滤掉 {record.id}，因为它缺少关键位点")

    ensure_core_sequence(core_fasta, filtered_fasta)
    print("✅ BLAST 结果检查完成，核心序列已保留")

# **5️⃣ 运行 CD-HIT，并检查关键位点，同时强制保留核心序列**
def run_cd_hit(input_fasta, output_fasta, key_sites, core_fasta):
    subprocess.run([
        "cd-hit", "-i", input_fasta, "-o", output_fasta, "-c", "0.95", "-n", "5", "-g", "1"
    ])

    missing_sites = check_validity(output_fasta, key_sites)
    if missing_sites:
        print(f"❌ CD-HIT 结果丢失关键位点 {sorted(missing_sites)}，终止执行！")
        exit(1)
    
    ensure_core_sequence(core_fasta, output_fasta)
    print("✅ CD-HIT 关键位点检查通过，核心序列已保留")

# **6️⃣ 运行 MUSCLE，并自动修复关键位点，同时强制保留核心序列**
def run_muscle(input_fasta, output_fasta, key_sites, core_fasta):
    try:
        result = subprocess.run(["muscle", "-version"], capture_output=True, text=True)
        if "5." in result.stdout:  
            subprocess.run(["muscle", "-align", input_fasta, "--output", output_fasta])
        else:  
            subprocess.run(["muscle", "-in", input_fasta, "-out", output_fasta])
        
        print("✅ MUSCLE 运行成功")

        missing_sites = check_validity(output_fasta, key_sites)
        if missing_sites:
            print(f"❌ MUSCLE 结果丢失关键位点 {sorted(missing_sites)}，尝试自动修复...")
            recover_muscle_alignment(output_fasta, key_sites)
            print("✅ 关键位点修复完成，已更新 MUSCLE 结果")
        
        ensure_core_sequence(core_fasta, output_fasta)
        print("✅ MUSCLE 关键位点检查通过，核心序列已保留")

    except Exception as e:
        print(f"❌ MUSCLE 运行失败: {e}")
        exit(1)
# **7️⃣ 关键位点检查函数**
def check_validity(fasta_file, key_sites):
    missing_sites = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        for pos in key_sites:
            if pos - 1 < len(seq) and seq[pos - 1] == "-":
                missing_sites.add(pos)

    return missing_sites

# **8️⃣ 自动修复关键位点**
def recover_muscle_alignment(msa_fasta, key_sites):
    records = list(SeqIO.parse(msa_fasta, "fasta"))

    for record in records:
        seq = list(str(record.seq))
        for pos in key_sites:
            if pos - 1 < len(seq) and seq[pos - 1] == "-":
                seq[pos - 1] = "X"  # **替换回 "X" 作为占位符**

        record.seq = Seq("".join(seq))

    with open(msa_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")

    print("✅ 关键位点修复完成，已更新 MUSCLE 结果")

# **9️⃣ 运行完整流程**
def main(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    json_file = os.path.join(input_dir, "BSH_t.json")
    core_fasta = os.path.join(input_dir, "BSH_t_core_it.fasta")
    homolog_fasta = os.path.join(input_dir, "BSH_homo.fasta")

    blast_output = os.path.join(output_dir, "blast_results.txt")
    filtered_fasta = os.path.join(output_dir, "BSH_t_filtered.fasta")
    cd_hit_output = os.path.join(output_dir, "BSH_t_cd_hit.fasta")
    muscle_output = os.path.join(output_dir, "BSH_t_msa.fasta")

    key_sites = parse_json(json_file)
    print(f"🔍 关键位点: {key_sites}")

    run_blast(core_fasta, homolog_fasta, blast_output, output_dir)
    print("✅ BLAST 完成")

    filter_blast_results(blast_output, homolog_fasta, key_sites, filtered_fasta, core_fasta)
    print("✅ 过滤 BLAST 结果，保留关键位点")

    run_cd_hit(filtered_fasta, cd_hit_output, key_sites, core_fasta)

    run_muscle(cd_hit_output, muscle_output, key_sites, core_fasta)

if __name__ == "__main__":
    input_dir = "/mnt/c/Users/ASUS/OneDrive/zby_2025/tsaas/nature/bai/Bai_hmm/BSH_t"
    output_dir = "/mnt/c/Users/ASUS/OneDrive/zby_2025/tsaas/nature/bai/Bai_hmm/BSH_t/output_results"
    main(input_dir, output_dir)
