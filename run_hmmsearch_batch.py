#conda activate hmm_env
#python run_hmmsearch_batch.py
import os
import subprocess

def run_prodigal(dna_path, protein_path):
    cmd = [
        "prodigal",
        "-i", dna_path,
        "-a", protein_path,
        "-p", "meta",  # 用于宏基因组或无注释基因组
        "-q"  # 静默运行
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"❌ Prodigal 失败: {dna_path}")
        print(result.stderr)
    else:
        print(f"✅ Prodigal 转换完成: {os.path.basename(protein_path)}")

def run_hmmsearch(hmm_model, protein_file, output_file):
    cmd = [
        "hmmsearch",
        "--tblout", output_file,
        "--noali",
        hmm_model,
        protein_file
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"❌ hmmsearch 失败: {protein_file}")
        print(result.stderr)
    else:
        print(f"✅ HMM 搜索完成: {os.path.basename(protein_file)}")

def main():
    hmm_model = "BaiE/BaiE_final.hmm"
    genome_dir = "test"
    protein_dir = os.path.join(genome_dir, "AA")
    output_dir = "BaiE/hmmsearch"
    os.makedirs(protein_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    for fname in os.listdir(genome_dir):
        if fname.endswith(".fna"):
            sample = os.path.splitext(fname)[0]
            dna_path = os.path.join(genome_dir, fname)
            protein_path = os.path.join(protein_dir, sample + ".faa")
            output_path = os.path.join(output_dir, sample + ".tbl")

            # Step 1: Prodigal 翻译 DNA → 蛋白
            run_prodigal(dna_path, protein_path)

            # Step 2: 使用蛋白序列进行 HMM 搜索
            run_hmmsearch(hmm_model, protein_path, output_path)

    print("🎉 所有序列翻译与搜索完成！")

if __name__ == "__main__":
    main()
