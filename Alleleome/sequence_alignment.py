# Alleleome generation step II
# Alignment of each allele sequnces (amino acid sequence) of each gene with consensus sequence using BLASTp.
import logging
import subprocess
from pathlib import Path
import gzip
import pandas as pd


def align_sequences(
    gene_list, out_dir, sequence_type="nucleotide", p=1
):
    logging.info("Starting alignment")
    out_dir = Path(out_dir)

    if p == 1:
        logging.info(f"Aligning sequences sequentially (p={p})")
        for gene_id in gene_list:
            align_single_gene(gene_id, out_dir, sequence_type=sequence_type)
    else:
        logging.info(f"Aligning sequences in parallel (p={p})")

        from multiprocessing import Pool
        from itertools import repeat

        chunksize = max(min(len(gene_list) // p, 500), 1)
        logging.info(f"Parallel chunksize = {chunksize}")
        with Pool(p) as pool:
            for _ in pool.imap_unordered(align_single_gene_parallel, zip(gene_list, repeat(out_dir), repeat(sequence_type)), chunksize=chunksize):
                pass
    
    logging.info("Completed  align_sequences in sequence_alignment")

def align_single_gene_parallel(args):
    gene_id, out_dir, sequence_type = args
    align_single_gene(gene_id, out_dir, sequence_type=sequence_type)

def align_single_gene(gene_id, out_dir, sequence_type="nucleotide"):
    out_dir = Path(out_dir)
    ext, blast = {"nucleotide": ("fna", "blastn"), "amino_acid": ("faa", "blastp")}[sequence_type]
    out_file = out_dir / "output" / gene_id / f"{sequence_type}_blast_out_{gene_id}.xml.gz"
    if out_file.is_file():
        logging.info(f"Outputs for {gene_id} already present, skipping.")
        return
    args = (
        blast,
        "-query",
        out_dir / "input" / gene_id / f"pan_genes.{ext}",
        "-subject",
        out_dir / "output" / gene_id / f"{sequence_type}_consensus_{gene_id}.{ext}",
        "-outfmt",
        "5",
    )
    with gzip.open(out_file, "wb") as f:
        subprocess.run(args, stdout=f)
    
