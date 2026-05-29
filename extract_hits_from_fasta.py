#conda activate hmm_env
#python extract_hits_from_fasta.py
import os
from Bio import SeqIO

def parse_tblout(tbl_file, evalue_threshold=1e-3):#Official manual recommends thresholds
    hits = set()
    with open(tbl_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            try:
                evalue = float(parts[4])
                if evalue <= evalue_threshold:
                    hits.add(parts[0])
            except ValueError:
                continue
    return hits

def extract_sequences(fasta_file, hit_ids, output_file):
    total = 0
    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in hit_ids:
                SeqIO.write(record, out_f, "fasta")
                total += 1
    print(f"✅ Extract {os.path.basename(output_file)}: {total} hits")

def main():
    tbl_dir = "BaiE/hmmsearch"
    genome_dir = "test/AA"
    output_dir = os.path.join(tbl_dir, "hits")
    os.makedirs(output_dir, exist_ok=True)

    for file in os.listdir(tbl_dir):
        if file.endswith(".tbl"):
            tbl_path = os.path.join(tbl_dir, file)
            genome_name = os.path.splitext(file)[0]
            genome_path = os.path.join(genome_dir, genome_name + ".faa")
            output_path = os.path.join(output_dir, genome_name + "_hits.fasta")

            if not os.path.exists(genome_path):
                print(f"⚠️ No corresponding genome file found: {genome_path}")
                continue

            hit_ids = parse_tblout(tbl_path)
            extract_sequences(genome_path, hit_ids, output_path)

    print("🎉 All matching sequences have been extracted！")

if __name__ == "__main__":
    main()
