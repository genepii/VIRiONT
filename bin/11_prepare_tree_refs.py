#!/usr/bin/env python3

"""

Sélectionne un sous-ensemble de séquences de référence pour l'arbre phylogénétique :

- Pour le VHD : inclut l'intégralité du panel de référence pour garantir une

  topologie complète et stable sans rupture de bootstrap.

- Pour le VHB (logique d'origine conservée) :

  * Ne garde que N sous-génotypes "squelette" par génotype (lettres A-I).

  * Si un sous-génotype d'un génotype donné est VALIDATED dans le run,

    TOUTES les références de ce génotype sont incluses.

"""

import sys

import os

import re

import argparse

import pandas as pd

from collections import defaultdict

from Bio import SeqIO



def parse_args():

    # Supporte à la fois les arguments nommés (--ref-fasta) et positionnels

    if len(sys.argv) > 1 and sys.argv[1].startswith('-'):

        parser = argparse.ArgumentParser(description="Sélection des références pour l'arbre phylogénétique")

        parser.add_argument("--ref-fasta", required=True, help="Fichier FASTA de référence")

        parser.add_argument("--summary-tsv", required=True, help="Fichier résumé clinique TSV")

        parser.add_argument("--out", dest="out_fasta", required=True, help="Fichier FASTA de sortie")

        parser.add_argument("--skeleton-n", dest="skeleton_n", type=int, default=2, help="Taille du squelette par défaut")

        args = parser.parse_args()

        return args.ref_fasta, args.summary_tsv, args.out_fasta, args.skeleton_n

    else:

        if len(sys.argv) < 4:

            print("Usage: 11_filter_tree_refs.py <ref_fasta> <summary_tsv> <out_fasta> [skeleton_n=2]")

            sys.exit(1)

        ref_fasta = sys.argv[1]

        summary_tsv = sys.argv[2]

        out_fasta = sys.argv[3]

        skeleton_n = int(sys.argv[4]) if len(sys.argv) > 4 else 2

        return ref_fasta, summary_tsv, out_fasta, skeleton_n



def main():

    ref_fasta, summary_tsv, out_fasta, skeleton_n = parse_args()



    records = list(SeqIO.parse(ref_fasta, 'fasta'))

    if not records:

        print(f"Erreur : Aucune séquence trouvée dans {ref_fasta}")

        sys.exit(1)



    ref_basename_upper = os.path.basename(ref_fasta).upper()

    is_vhd = ("HDV" in ref_basename_upper) or ("VHD" in ref_basename_upper)



    # =========================================================================

    # BRANCHE VHD : CONSERVATION DE TOUTES LES RÉFÉRENCES DU PANEL

    # =========================================================================

    if is_vhd:

        print(f"[VHD détecté] Conservation de l'ensemble du panel de référence ({len(records)} séquences).")

        SeqIO.write(records, out_fasta, 'fasta')

        print(f" {len(records)} séquences de référence écrites dans {out_fasta} pour le VHD.")

        return



    # =========================================================================

    # BRANCHE VHB : LOGIQUE HISTORIQUE STRICTE (SQUELETTE & ZOOM RUN)

    # =========================================================================

    pattern = re.compile(r'^([A-I])(\d*)')



    # 1. Génotypes réellement validés dans ce run

    present_letters = set()

    try:

        if os.path.exists(summary_tsv) and os.path.getsize(summary_tsv) > 0:

            # Gestion séparateur TSV ou CSV

            with open(summary_tsv, 'r', encoding='utf-8') as f:

                sep = '\t' if '\t' in f.readline() else ';'

            df = pd.read_csv(summary_tsv, sep=sep)

            

            status_col = next((c for c in df.columns if c.lower() == 'status'), None)

            geno_col   = next((c for c in df.columns if c.lower() == 'genotype'), None)



            if status_col and geno_col:

                validated = df.loc[df[status_col].astype(str).str.upper() == 'VALIDATED', geno_col].astype(str)

                for g in validated:

                    m = pattern.match(g.strip())

                    if m:

                        present_letters.add(m.group(1))

    except Exception as e:

        print(f"Erreur lecture {summary_tsv}: {e}")



    print(f"[VHB] Génotypes validés détectés dans ce run : {sorted(present_letters) if present_letters else 'aucun'}")



    # 2. Regroupement des références par génotype (lettre A-I) et sous-génotype

    groups = defaultdict(lambda: defaultdict(list))

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



    # 3. Sélection : squelette (N sous-génotypes) par défaut, tout si le génotype est présent

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

            selected.append(subgeno_dict[sg][0])



    SeqIO.write(selected, out_fasta, 'fasta')

    print(f"\n [VHB] {len(selected)} séquences de référence sélectionnées pour l'arbre (sur {len(records)} disponibles).")



if __name__ == "__main__":

    main()

