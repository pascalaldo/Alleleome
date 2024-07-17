import logging
from pathlib import Path
import shutil

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_genbank_files(df_gene_presence_locustag, gbk_folder):
    """
    Parse GenBank files.

    Parameters:
    df_gene_presence_locustag (DataFrame): DataFrame with gene presence information
    gbk_folder (Path): Path to folder containing GenBank files

    Returns:
    DataFrame: DataFrame with parsed GenBank data
    """

    all_locustag_list = []
    for genome_id in list(df_gene_presence_locustag.columns):
        genbank_file_path = Path(gbk_folder) / f"{genome_id}.gbk"
        for record in SeqIO.parse(genbank_file_path, "genbank"):
            for feature in record.features:
                genome_data_list = []
                tag = feature.qualifiers.get("locus_tag")
                if not tag:
                    continue
                gene_id = df_gene_presence_locustag.index[df_gene_presence_locustag[genome_id] == tag[0]]
                if len(gene_id) != 1:
                    continue
                genome_data_list.append(tag[0])  # Locus tag
                genome_data_list.append(genome_id)  # Genome ID
                genome_data_list.append(gene_id[0])
                
                if "product" in feature.qualifiers.keys():
                    genome_data_list.append(
                        feature.qualifiers["product"][0]
                    )  # Prokka annotation
                else:
                    genome_data_list.append("")
                genome_data_list.append(
                    int(feature.location.start)
                )  # Start position
                genome_data_list.append(int(feature.location.end))  # End position

                # nuc_record = SeqRecord(
                #     feature.extract(record.seq),
                #     id=genome_data_list[0], # locus tag
                #     description=f"{(genome_data_list[3] if genome_data_list[3] != '' else 'Unknown')} | {genome_data_list[1]}",
                # )
                # with open(tmp_folder / f"{counter}.fna", "w") as f:
                #     SeqIO.write(nuc_record, f, "fasta")

                # genome_data_list.append(len(nuc_record))
                genome_data_list.append(len(feature))


                # aa_record = SeqRecord(
                #     feature.qualifiers["translation"][0] if "translation" in feature.qualifiers.keys() else "",
                #     id=genome_data_list[0], # locus tag
                #     description=f"{(genome_data_list[3] if genome_data_list[3] != '' else 'Unknown')} | {genome_data_list[1]}",
                # )
                # with open(tmp_folder / f"{counter}.faa", "w") as f:
                #     SeqIO.write(aa_record, f, "fasta")

                all_locustag_list.append(genome_data_list)

    all_locustag_df = pd.DataFrame(
        all_locustag_list,
        columns=[
            "Locus_Tag",
            "Genome_ID",
            "Gene_ID",
            "Prokka_Annotation",
            "Start_Position",
            "End_Position",
            "Nucleotide_Len",
        ],
    )
    # Use genome + locus tag as index, since there may be duplicate locus tags accross genomes.
    all_locustag_df.index = all_locustag_df["Genome_ID"] + "@" + all_locustag_df["Locus_Tag"] 
    return all_locustag_df

# def process_gene(
#     gene_id, sel_locustag_df, fna_path, faa_path, gbk_folder
# ):
#     logging.info(f"   Processing gene: {gene_id}")

#     filtered_df = sel_locustag_df[sel_locustag_df["Gene_ID"] == gene_id]

#     files = [(tmp_folder / f"{int(file_id)}") for file_id in filtered_df["File_ID"]]

#     for seq_type, out_file in [(".fna", fna_path), (".faa", faa_path)]:
#         with open(out_file, "w") as f_out:
#             for in_file_base in files:
#                 in_file = in_file_base.with_suffix(seq_type)
#                 with open(in_file, "r") as f_in:
#                     s = f_in.readlines()
#                     f_out.writelines(s)
#     logging.info(f"   Finished gene: {gene_id}")

def process_selected_genes(all_locustag_df, locustag_list, gene_list, out_dir, gbk_folder):
    out_dir = Path(out_dir)
    tmp_folder = Path(tmp_folder)
    sel_locustag_df = all_locustag_df.loc[locustag_list]
    
    for gene in gene_list:
        gene_dir = out_dir / "input" / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        fna_path = gene_dir / "pan_genes.fna"
        faa_path = gene_dir / "pan_genes.faa"
        if fna_path.is_file():
            fna_path.unlink()
        if faa_path.is_file():
            faa_path.unlink()
        # process_gene(gene, sel_locustag_df, fna_path, faa_path, tmp_folder)
    genome_ids = sel_locustag_df["Genome_ID"].unique()
    for i, genome_id in enumerate(genome_ids):
        print(f"Writing genome #{i+1}/{len(genome_ids)} to fasta.")
        genbank_file_path = Path(gbk_folder) / f"{genome_id}.gbk"
        for record in SeqIO.parse(genbank_file_path, "genbank"):
            for feature in record.features:
                tag = feature.qualifiers.get("locus_tag")
                if not tag:
                    continue
                locustag = tag[0]
                try:
                    locustag_info = sel_locustag_df.loc[f"{genome_id}@{locustag}", :]
                except:
                    continue
                gene_id = str(locustag_info["Gene_ID"])
                if not gene_id in gene_list:
                    continue
                
                annotation = locustag_info["Prokka_Annotation"]
                if annotation == "":
                    annotation = "Unknown"

                gene_dir = out_dir / "input" / gene_id
                nuc_record = SeqRecord(
                    feature.extract(record.seq),
                    id=locustag, # locus tag
                    description=f"{annotation} | {genome_id}",
                )
                with open(gene_dir / f"pan_genes.fna", "a") as f:
                    SeqIO.write(nuc_record, f, "fasta")
                
                aa_record = SeqRecord(
                    feature.qualifiers["translation"][0] if "translation" in feature.qualifiers.keys() else "",
                    id=locustag, # locus tag
                    description=f"{annotation} | {genome_id}",
                )
                with open(gene_dir / f"pan_genes.faa", "a") as f:
                    SeqIO.write(aa_record, f, "fasta")
