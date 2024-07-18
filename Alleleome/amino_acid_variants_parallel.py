# Alleleome generation step III- Parsing the results and generating the amino acid mutations
import logging
from . import amino_acid_variants_parallel_helper
from itertools import repeat
from multiprocessing import Pool
import pandas as pd


def generate_amino_acid_vars(gene_list, out_dir, aa_vars_path, p=8):
    gene_list_len = len(gene_list)
    chunksize = max(min(gene_list_len // p, 50), 1)
    logging.info(f"Parallel chunksize = {chunksize}")
    counter = 0

    with open(aa_vars_path, "w") as f:
        with Pool(p) as pool:
            for result in pool.imap_unordered(
                amino_acid_variants_parallel_helper.generate_amino_acid_vars,
                zip(gene_list, repeat(out_dir)),
                chunksize=chunksize,
            ):
                logging.info(f"Processing AAV result of gene #{counter+1}/{gene_list_len}")
                df = pd.DataFrame(result)
                df.to_csv(f, header=(counter == 0), index=False)
                counter += 1
