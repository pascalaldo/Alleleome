# Alleleome generation step II
# Alignment of each allele sequnces (amino acid sequence) of each gene with consensus sequence using BLASTp.
import logging
import subprocess
from pathlib import Path

import pandas as pd


def align_sequences(
    gene_list, out_dir, sequence_type="nucleotide"
):
    logging.info("Starting alignment")
    out_dir = Path(out_dir)

    for gene_id in gene_list:
        align_single_gene(gene_id, out_dir / gene_id, sequence_type=sequence_type)
    
    logging.info("Completed  align_sequences in sequence_alignment")

def align_single_gene(gene_id, out_dir, sequence_type="nucleotide"):
    out_dir = Path(out_dir)
    ext, blast = {"nucleotide": ("fna", "blastn"), "amino_acid": ("faa", "blastp")}[sequence_type]
    out_file_name = out_dir / "output" / f"{sequence_type}_blast_out_{gene_id}.xml"
    args = (
        blast,
        "-query",
        out_dir / "input" / f"pan_genes.{ext}",
        "-subject",
        out_dir / "output" / f"{sequence_type}_consensus_{gene_id}.{ext}",
        "-outfmt",
        "5",
    )
    with open(out_file_name, "w+") as outfile:
        subprocess.run(args, stdout=outfile)
    
