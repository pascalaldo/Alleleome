import logging
from pathlib import Path

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def process_gene(
    gene_id, sel_locustag_df, fna_path, faa_path
):
    logging.info(f"   Processing gene: {gene_id}")

    filtered_df = sel_locustag_df[sel_locustag_df["Gene_ID"] == gene_id]
    nucleotide_records = [
        SeqRecord(
            Seq(row["Nucleotide_Seq"]),
            id=row["Locus_Tag"],
            description=f"{row.get('Prokka_Annotation', 'Unknown')} | {row['Genome_ID']}",
        )
        for idx, row in filtered_df.iterrows()
    ]
    amino_acid_records = [
        SeqRecord(
            Seq(row["Amino_Acid_Seq"]),
            id=row["Locus_Tag"],
            description=f"{row.get('Prokka_Annotation', 'Unknown')} | {row['Genome_ID']}",
        )
        for idx, row in filtered_df.iterrows()
    ]
    with open(fna_path, "w") as fasta_n_file:
        SeqIO.write(nucleotide_records, fasta_n_file, "fasta")
    with open(faa_path, "w") as fasta_aa_file:
        SeqIO.write(amino_acid_records, fasta_aa_file, "fasta")
    logging.info(f"   Finished gene: {gene_id}")

def process_selected_genes(all_locustag_df, locustag_list, gene_list, out_dir):
    out_dir = Path(out_dir)
    sel_locustag_df = all_locustag_df.loc[locustag_list]
    
    for gene in gene_list:
        gene_dir = out_dir / "input" / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        fna_path = gene_dir / "pan_genes.fna"
        faa_path = gene_dir / "pan_genes.faa"
        if fna_path.is_file() and faa_path.is_file():
            continue
        process_gene(gene, sel_locustag_df, fna_path, faa_path)
