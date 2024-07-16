import logging
from pathlib import Path
import shutil

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def process_gene(
    gene_id, sel_locustag_df, fna_path, faa_path, tmp_folder
):
    logging.info(f"   Processing gene: {gene_id}")

    filtered_df = sel_locustag_df[sel_locustag_df["Gene_ID"] == gene_id]

    files = [(tmp_folder / f"{int(file_id)}") for file_id in filtered_df["File_ID"]]

    for seq_type, out_file in [(".fna", fna_path), (".faa", faa_path)]:
        with open(out_file, "w") as f_out:
            for in_file_base in files:
                in_file = in_file_base.with_suffix(seq_type)
                with open(in_file, "r") as f_in:
                    s = f_in.readlines()
                    f_out.writelines(s)
    logging.info(f"   Finished gene: {gene_id}")

def process_selected_genes(all_locustag_df, locustag_list, gene_list, out_dir, tmp_folder):
    out_dir = Path(out_dir)
    tmp_folder = Path(tmp_folder)
    sel_locustag_df = all_locustag_df.loc[locustag_list]
    
    for gene in gene_list:
        gene_dir = out_dir / "input" / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        fna_path = gene_dir / "pan_genes.fna"
        faa_path = gene_dir / "pan_genes.faa"
        if fna_path.is_file() and faa_path.is_file():
            continue
        process_gene(gene, sel_locustag_df, fna_path, faa_path, tmp_folder)
