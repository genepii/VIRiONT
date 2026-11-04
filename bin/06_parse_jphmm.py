#!/usr/bin/env python3
import sys
import os
import re

NON_GENOTYPE_TOKENS = {
    "3'-INSERTION", "5'-INSERTION", "INSERTION", "UR", 
    "UNASSIGNED", "UNKNOWN", "GAP", "DELETION"
}

def parse_jphmm(output_dir, sample_id):
    summary_tsv = f"{sample_id}_recombination_summary.tsv"
    
    target_file = None
    for root, dirs, files in os.walk(output_dir):
        for f in files:
            if f == "recombination.txt" or f.endswith("_recombination.txt"):
                target_file = os.path.join(root, f)
                break
        if target_file:
            break

    segments = []
    viral_subtypes = []

    if target_file and os.path.exists(target_file):
        with open(target_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith(">"):
                    continue
                parts = re.split(r"\s+", line)
                if len(parts) >= 3 and parts[0].isdigit() and parts[1].isdigit():
                    start, end, raw_st = parts[0], parts[1], parts[2]
                    segments.append(f"{raw_st}:{start}-{end}")
                    
                    clean_st = raw_st.strip()
                    if clean_st.upper() not in NON_GENOTYPE_TOKENS:
                        geno_clade = clean_st[0].upper()
                        if geno_clade in "ABCDEFGHI" and clean_st not in viral_subtypes:
                            viral_subtypes.append(clean_st)

    distinct_clades = {st[0].upper() for st in viral_subtypes if st[0].upper() in "ABCDEFGHI"}
    is_recombinant = "YES" if len(distinct_clades) > 1 else "NO"
    
    subtypes_str = "+".join(viral_subtypes) if viral_subtypes else "Undetermined"
    breakpoints_str = ";".join(segments) if segments else "None"

    with open(summary_tsv, "w") as out:
        out.write("sample_id\tis_recombinant\tsubtypes\tbreakpoints\n")
        out.write(f"{sample_id}\t{is_recombinant}\t{subtypes_str}\t{breakpoints_str}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("Usage: 06_parse_jphmm.py <jphmm_output_dir> <sample_id>")
    parse_jphmm(sys.argv[1], sys.argv[2])