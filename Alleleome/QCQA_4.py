import logging
import os
import shutil
from pathlib import Path

import pandas as pd


def process_genes(
    pangenome_alignments_dir_path, alleleome_dir_path, pangene_summary_csv=None, pan_core="Core"
):
    """
    Processes core genes by copying gene sequence files and creating output directories for each gene.

    Parameters:
    - main_dir: Path to the main directory containing gene folders.
    - alleleome_path: Path to the directory containing the gene data CSV files.
    """
    try:
        logging.info("Starting process_genes in QCQA_4")
        pangenome_alignments_dir_path = Path(pangenome_alignments_dir_path)
        # Read genes data
        alleleome_dir_path = Path(alleleome_dir_path)
        alleleome_dir_path.mkdir(parents=True, exist_ok=True)

        if pan_core == "Core":
            if pangene_summary_csv is None:
                pangene_summary_csv = alleleome_dir_path / "df_pangene_summary_v2.csv"
            else:
                pangene_summary_csv = Path(pangene_summary_csv)

            assert (
                pangene_summary_csv.is_file()
            ), f"Cannot find pangene_summary table at {pangene_summary_csv}"

            df = pd.read_csv(pangene_summary_csv)
            gene_list = df[df["pangenome_class_2"] == "Core"]["Gene"].unique().tolist()
        elif pan_core == "Pan":
            # Read genes data
            df = pd.read_csv(os.path.join(alleleome_dir_path,'nuc_genes_present_in_above_5_percent_of_strains.csv'))
            gene_list = df['Gene'].tolist()
        else:
            raise ValueError("Unrecognized alleleome type, should be Core or Pan.")

        # Read core gene alleles data
        df1 = pd.read_csv(
            os.path.join(
                alleleome_dir_path,
                "alleles_with_length_less_than_2std_less_than_mean_length.csv",
            )
        )
        edit_list = df1["Gene"].tolist()
        new_gene_list = set(edit_list)

        # Process each core gene
        for gene in gene_list:
            # Copy gene sequence files if gene not in new_gene_list
            if gene not in new_gene_list:
                aa_allele_path = pangenome_alignments_dir_path / gene / "input"
                aa_file = aa_allele_path / "pangenes.faa"
                new_file = aa_allele_path / "pan_genes.faa"
                na_file = aa_allele_path / "pangenes.fna"
                new_na_file = aa_allele_path / "pan_genes.fna"
                shutil.copyfile(aa_file, new_file)
                shutil.copyfile(na_file, new_na_file)

            # Create output directory for the gene
            output_path = pangenome_alignments_dir_path / gene / "output"
            output_path.mkdir(exist_ok=True, parents=True)
        logging.info("Completed process_genes in QCQA_4")
    except Exception as e:
        logging.error(f"Error in process_genes in QCQA_4: {e}")
        raise
