# conda activate cdhit_env
# python key_point.py
# cat BaiN_ncbi.fasta BaiN_uniprot.fasta > BaiN_homo.fasta
import json
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# **1️⃣ Parse *.json files and extract key sites**
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

# **2️⃣ Ensure the core sequence *_core always exists.**
def ensure_core_sequence(core_fasta, output_fasta):
    core_seq_record = list(SeqIO.parse(core_fasta, "fasta"))[0]  # read the core sequence

    # Read all sequences in the current file
    records = list(SeqIO.parse(output_fasta, "fasta"))
    record_ids = [record.id for record in records]

    # If the core sequence is not in the results, add it.
    if core_seq_record.id not in record_ids:
        print(f"⚠ Forced preservation of core sequences {core_seq_record.id}")
        records.append(core_seq_record)

    # Rewrite the file to ensure the core sequence is preserved.
    with open(output_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")

# **3️⃣ Run BLAST**
def run_blast(query_fasta, db_fasta, output_file, output_dir):
    db_name = os.path.join(output_dir, "blast_db")
    
    subprocess.run(["makeblastdb", "-in", db_fasta, "-dbtype", "prot", "-out", db_name])
    
    subprocess.run([
        "blastp", "-query", query_fasta, "-db", db_name, 
        "-out", output_file, "-evalue", "1e-5", "-max_target_seqs", "500", "-outfmt", "6"
    ])

# **4️⃣ Filter the BLAST results and check key sites, while forcibly preserving the core sequence**
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
                    print(f"⚠ Filter out {record.id}，because it lacks key sites")

    ensure_core_sequence(core_fasta, filtered_fasta)
    print("✅ BLAST result check complete, core sequence preserved.")

# **5️⃣ Run CD-HIT and check key sites while forcibly preserving the core sequence.**
def run_cd_hit(input_fasta, output_fasta, key_sites, core_fasta):
    subprocess.run([
        "cd-hit", "-i", input_fasta, "-o", output_fasta, "-c", "0.95", "-n", "5", "-g", "1"
    ])

    missing_sites = check_validity(output_fasta, key_sites)
    if missing_sites:
        print(f"❌ CD-HIT results showed loss of key sites. {sorted(missing_sites)}，stop！")
        exit(1)
    
    ensure_core_sequence(core_fasta, output_fasta)
    print("✅ CD-HIT key site check passed, core sequence preserved.")

# **6️⃣ Run MUSCLE and automatically repair critical sites while forcibly preserving the core sequence**
def run_muscle(input_fasta, output_fasta, key_sites, core_fasta):
    try:
        result = subprocess.run(["muscle", "-version"], capture_output=True, text=True)
        if "5." in result.stdout:  
            subprocess.run(["muscle", "-align", input_fasta, "--output", output_fasta])
        else:  
            subprocess.run(["muscle", "-in", input_fasta, "-out", output_fasta])
        
        print("✅ MUSCLE ran successfully")

        missing_sites = check_validity(output_fasta, key_sites)
        if missing_sites:
            print(f"❌ MUSCLE results show loss of key sites {sorted(missing_sites)}，attempt automatic repair...")
            recover_muscle_alignment(output_fasta, key_sites)
            print("✅ Key site repair completed, MUSCLE results updated")
        
        ensure_core_sequence(core_fasta, output_fasta)
        print("✅ The MUSCLE key site check passed, and the core sequence has been preserved")

    except Exception as e:
        print(f"❌ MUSCLE failed: {e}")
        exit(1)
# **7️⃣ Key site checking function**
def check_validity(fasta_file, key_sites):
    missing_sites = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        for pos in key_sites:
            if pos - 1 < len(seq) and seq[pos - 1] == "-":
                missing_sites.add(pos)

    return missing_sites

# **8️⃣ Automatic repair of critical sites**
def recover_muscle_alignment(msa_fasta, key_sites):
    records = list(SeqIO.parse(msa_fasta, "fasta"))

    for record in records:
        seq = list(str(record.seq))
        for pos in key_sites:
            if pos - 1 < len(seq) and seq[pos - 1] == "-":
                seq[pos - 1] = "X"  # **Replace "X" as a placeholder**

        record.seq = Seq("".join(seq))

    with open(msa_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")

    print("✅ Key site repair completed, MUSCLE results updated")

# **9️⃣ Run the complete process**
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
    print(f"🔍 key sites: {key_sites}")

    run_blast(core_fasta, homolog_fasta, blast_output, output_dir)
    print("✅ BLAST done")

    filter_blast_results(blast_output, homolog_fasta, key_sites, filtered_fasta, core_fasta)
    print("✅ Filter BLAST results and retain key sites")

    run_cd_hit(filtered_fasta, cd_hit_output, key_sites, core_fasta)

    run_muscle(cd_hit_output, muscle_output, key_sites, core_fasta)

if __name__ == "__main__":
    input_dir = "/mnt/c/Users/ASUS/OneDrive/zby_2025/tsaas/nature/bai/Bai_hmm/BSH_t"
    output_dir = "/mnt/c/Users/ASUS/OneDrive/zby_2025/tsaas/nature/bai/Bai_hmm/BSH_t/output_results"
    main(input_dir, output_dir)
