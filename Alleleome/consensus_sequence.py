# Alleleome generation step I- Get the consensus of the allels (amino acid and nucleotide sequences) of each core gene
import logging
import os
from pathlib import Path

import pandas as pd
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MafftCommandline


def build_consensus(
    gene_list, out_dir
):
    logging.info("Starting build_consensus in build_consensus_sequence")
    out_dir = Path(out_dir)

    for gene_id in gene_list:
        build_single_gene_consensus(gene_id, out_dir / gene_id)
    
    logging.info("Completed build_consensus in build_consensus_sequence")

def build_single_gene_consensus(gene_id, out_dir):
    out_dir = Path(out_dir)
    alignment_input_dir = out_dir / "input"
    alignment_output_dir = out_dir / "output"
    alignment_output_dir.mkdir(parents=True, exist_ok=True)

    for cons_type, ext in {"amino_acid": "faa", "nucleotide": "fna"}.items():
        generate_consensus(
            gene_id,
            alignment_input_dir / f"pan_genes.{ext}",
            alignment_output_dir / f"mafft_{cons_type}_{gene_id}.fasta",
            alignment_output_dir / f"{cons_type}_consensus_{gene_id}.{ext}"
        )


def generate_consensus(gene_id, input_path, mafft_output_path, consensus_output_path):
    if os.stat(input_path).st_size == 0:
        print("The file is empty")
        return

    mafft_cline = MafftCommandline(input=input_path)
    stdout, stderr = mafft_cline()

    with open(mafft_output_path, "w") as handle:
        handle.write(stdout)

    myalign = AlignIO.read(mafft_output_path, "fasta")
    summary = AlignInfo.SummaryInfo(myalign)
    consensus = summary.dumb_consensus(threshold=0.5)
    seq = str(consensus).upper()

    with open(consensus_output_path, "w") as consensus_file:
        consensus_seq = "".join([">" + gene_id + "\n" + seq])
        consensus_file.write(consensus_seq)