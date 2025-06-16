# conda activate cdhit_env
# python no_key_point.py
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# **1️⃣ 确保核心序列 *_core 始终存在**
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

# **2️⃣ 运行 BLAST**
def run_blast(query_fasta, db_fasta, output_file, output_dir):
    db_name = os.path.join(output_dir, "blast_db")
    
    subprocess.run(["makeblastdb", "-in", db_fasta, "-dbtype", "prot", "-out", db_name])
    
    subprocess.run([
        "blastp", "-query", query_fasta, "-db", db_name, 
        "-out", output_file, "-evalue", "1e-5", "-max_target_seqs", "500", "-outfmt", "6"
    ])

# **3️⃣ 过滤 BLAST 结果，同时强制保留核心序列**
def filter_blast_results(blast_output, db_fasta, filtered_fasta, core_fasta):
    valid_ids = set()
    
    with open(blast_output, "r") as f:
        for line in f:
            parts = line.split("\t")
            seq_id = parts[1]
            valid_ids.add(seq_id)

    with open(filtered_fasta, "w") as f_out:
        for record in SeqIO.parse(db_fasta, "fasta"):
            if record.id in valid_ids:
                SeqIO.write(record, f_out, "fasta")

    ensure_core_sequence(core_fasta, filtered_fasta)
    print("✅ BLAST 结果检查完成，核心序列已保留")

# **4️⃣ 运行 CD-HIT，同时强制保留核心序列**
def run_cd_hit(input_fasta, output_fasta, core_fasta):
    subprocess.run([
        "cd-hit", "-i", input_fasta, "-o", output_fasta, "-c", "0.95", "-n", "5", "-g", "1"
    ])
    
    ensure_core_sequence(core_fasta, output_fasta)
    print("✅ CD-HIT 运行完成，核心序列已保留")

# **5️⃣ 运行 MUSCLE，同时强制保留核心序列**
def run_muscle(input_fasta, output_fasta, core_fasta):
    try:
        result = subprocess.run(["muscle", "-version"], capture_output=True, text=True)
        if "5." in result.stdout:  
            subprocess.run(["muscle", "-align", input_fasta, "--output", output_fasta])
        else:  
            subprocess.run(["muscle", "-in", input_fasta, "-out", output_fasta])
        
        print("✅ MUSCLE 运行成功")
        
        ensure_core_sequence(core_fasta, output_fasta)
        print("✅ MUSCLE 关键位点检查通过，核心序列已保留")

    except Exception as e:
        print(f"❌ MUSCLE 运行失败: {e}")
        exit(1)

# **6️⃣ 运行完整流程**
def main(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    core_fasta = os.path.join(input_dir, "BaiN_core.fasta")
    homolog_fasta = os.path.join(input_dir, "BaiN_homo.fasta")

    blast_output = os.path.join(output_dir, "blast_results.txt")
    filtered_fasta = os.path.join(output_dir, "BaiN_filtered.fasta")
    cd_hit_output = os.path.join(output_dir, "BaiN_cd_hit.fasta")
    muscle_output = os.path.join(output_dir, "BaiN_msa.fasta")

    run_blast(core_fasta, homolog_fasta, blast_output, output_dir)
    print("✅ BLAST 完成")

    filter_blast_results(blast_output, homolog_fasta, filtered_fasta, core_fasta)

    run_cd_hit(filtered_fasta, cd_hit_output, core_fasta)

    run_muscle(cd_hit_output, muscle_output, core_fasta)

if __name__ == "__main__":
    input_dir = "/mnt/c/Users/ASUS/OneDrive/zby_2025/tsaas/nature/bai/Bai_hmm/BaiN"
    output_dir = "/mnt/c/Users/ASUS/OneDrive/zby_2025/tsaas/nature/bai/Bai_hmm/BaiN/output_results"
    main(input_dir, output_dir)
