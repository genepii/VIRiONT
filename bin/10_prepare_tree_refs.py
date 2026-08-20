#!/usr/bin/env python3
"""
Sélectionne un sous-ensemble de séquences de référence pour l'arbre phylogénétique :
- Par défaut, ne garde que N sous-génotypes "squelette" par génotype (letter A-I),
  pour aérer l'arbre.
- Si un sous-génotype d'un génotype donné est VALIDATED dans le run (présent dans
  SUMMARY_Multi_Infection.tsv), TOUTES les références de ce génotype sont incluses,
  pour donner la meilleure résolution phylogénétique autour des variants détectés.
"""
import sys
import re
import pandas as pd
from collections import defaultdict
from Bio import SeqIO

def main():
    if len(sys.argv) < 4:
        print("Usage: 10_prepare_tree_refs.py <ref_fasta> <summary_tsv> <out_fasta> [skeleton_n=2]")
        sys.exit(1)

    ref_fasta   = sys.argv[1]
    summary_tsv = sys.argv[2]
    out_fasta   = sys.argv[3]
    skeleton_n  = int(sys.argv[4]) if len(sys.argv) > 4 else 2

    pattern = re.compile(r'^([A-I])(\d*)')

    # 1. Génotypes réellement validés dans ce run
    present_letters = set()
    try:
        df = pd.read_csv(summary_tsv, sep='\t')
        if 'status' in df.columns and 'genotype' in df.columns:
            validated = df.loc[df['status'] == 'VALIDATED', 'genotype'].astype(str)
            for g in validated:
                m = pattern.match(g.strip())
                if m:
                    present_letters.add(m.group(1))
    except Exception as e:
        print(f"Erreur lecture {summary_tsv}: {e}")

    print(f"Génotypes validés détectés dans ce run : {sorted(present_letters) if present_letters else 'aucun'}")

    # 2. Regroupement des références par génotype (letter) et sous-génotype
    records = list(SeqIO.parse(ref_fasta, 'fasta'))
    groups = defaultdict(lambda: defaultdict(list))  # letter -> subgeno -> [records]
    for rec in records:
        token = rec.id.split()[0]
        m = pattern.match(token)
        if m:
            letter = m.group(1)
            num = m.group(2)
            subgeno = letter + num if num else letter
        else:
            letter = 'OTHER'
            subgeno = token
        groups[letter][subgeno].append(rec)

    # 3. Sélection : squelette (N sous-génotypes) par défaut, tout si le génotype
    #    est présent dans le run
    selected = []
    for letter, subgeno_dict in groups.items():
        subgeno_list = sorted(subgeno_dict.keys())
        if letter in present_letters:
            chosen = subgeno_list
            print(f"  Génotype {letter} : PRÉSENT dans le run -> {len(chosen)} sous-génotype(s) inclus ({', '.join(chosen)})")
        else:
            chosen = subgeno_list[:skeleton_n]
            print(f"  Génotype {letter} : absent du run -> squelette de {len(chosen)} sous-génotype(s) ({', '.join(chosen)})")
        for sg in chosen:
            # une seule séquence représentative par sous-génotype (évite les doublons
            # si plusieurs souches de référence existent pour le même code)
            selected.append(subgeno_dict[sg][0])

    SeqIO.write(selected, out_fasta, 'fasta')
    print(f"\n {len(selected)} séquences de référence sélectionnées pour l'arbre (sur {len(records)} disponibles).")


if __name__ == "__main__":
    main()