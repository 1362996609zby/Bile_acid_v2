import os
import json
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#conda activate hmm_env
#python build_annotated_hmm.py
def load_fasta(fasta_file):
    return {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}

def parse_annotation_file(json_file):
    with open(json_file) as f:
        data = json.load(f)
    
    annotations = {}

    if "features" in data:  # BaiA1格式
        seq_id = data.get("primaryAccession", "unknown")
        annotations[seq_id] = {"AS": [], "BS": []}
        for feature in data["features"]:
            label = "AS" if "Active" in feature["type"] else "BS"
            start = feature["location"]["start"]["value"]
            end = feature["location"]["end"]["value"]
            annotations[seq_id][label].extend(range(start, end + 1))

    else:  # BaiB结构域预测样式
        for domain_id, seq_map in data.items():
            for seq_id, marks in seq_map.items():
                if seq_id not in annotations:
                    annotations[seq_id] = {"AS": [], "BS": []}
                annotations[seq_id]["AS"].extend(marks.get("AS", []))
                annotations[seq_id]["BS"].extend(marks.get("BS", []))

    return annotations

def write_sto_with_annotations(fasta_dict, annotation_dict, output_sto):
    with open(output_sto, "w") as out:
        out.write("# STOCKHOLM 1.0\n")
        for seq_id, record in fasta_dict.items():
            seq = str(record.seq)
            out.write(f"{seq_id} {seq}\n")

            for label in ["AS", "BS"]:
                annot = ["."] * len(seq)
                for pos in annotation_dict.get(seq_id, {}).get(label, []):
                    if 1 <= pos <= len(seq):
                        annot[pos - 1] = "*" if label == "AS" else "+"
                out.write(f"#=GR {seq_id} {label} {''.join(annot)}\n")
        out.write("//\n")
    print(f"✅ 已生成注释对齐文件: {output_sto}")

def build_hmm(sto_file, hmm_file):
    cmd = ["hmmbuild", "--amino", hmm_file, sto_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"❌ hmmbuild failed:\n{result.stderr}")
    print(f"✅ 已成功构建 HMM 模型: {hmm_file}")

def main():
    # === 用户可自定义 ===
    fasta_file = "BSH_t/BSH_t_msa.fasta"
    json_files = ["BSH_t/BSH_t.json"]
    sto_output = "BSH_t/BSH_t_final.sto"
    hmm_output = "BSH_t/BSH_t_final.hmm"

    # === 加载序列和注释 ===
    seqs = load_fasta(fasta_file)
    all_annots = {}

    for jf in json_files:
        ann = parse_annotation_file(jf)
        for sid in ann:
            if sid not in all_annots:
                all_annots[sid] = {"AS": [], "BS": []}
            all_annots[sid]["AS"].extend(ann[sid]["AS"])
            all_annots[sid]["BS"].extend(ann[sid]["BS"])

    # 去重
    for sid in all_annots:
        all_annots[sid]["AS"] = sorted(set(all_annots[sid]["AS"]))
        all_annots[sid]["BS"] = sorted(set(all_annots[sid]["BS"]))
    # ✅ 若没有注释文件，跳过注释加载
    if json_files:
        for jf in json_files:
            ann = parse_annotation_file(jf)
            for sid in ann:
                if sid not in all_annots:
                    all_annots[sid] = {"AS": [], "BS": []}
                all_annots[sid]["AS"].extend(ann[sid]["AS"])
                all_annots[sid]["BS"].extend(ann[sid]["BS"])

    # 去重
        for sid in all_annots:
            all_annots[sid]["AS"] = sorted(set(all_annots[sid]["AS"]))
            all_annots[sid]["BS"] = sorted(set(all_annots[sid]["BS"]))
    # === 生成对齐注释并构建 HMM ===
    write_sto_with_annotations(seqs, all_annots, sto_output)
    build_hmm(sto_output, hmm_output)

if __name__ == "__main__":
    main()
