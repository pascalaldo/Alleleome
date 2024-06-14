import pandas as pd
from Bio import AlignIO
from pathlib import Path
from collections import Counter
import logging

def find_dominant_aa(gene_list, out_dir, dominant_aa_path):
    logging.info("Starting: preplot: find_dominant_aa")
    out_dir = Path(out_dir)
    # Initialize a list to hold consensus data
    consensus_data = []

    # Process each alignment file
    for gene in gene_list:
        msa_file_name = out_dir / "output" / gene / f'mafft_amino_acid_{gene}.fasta'
        # logging.info(msa_file_name)
        aln = AlignIO.read(msa_file_name, 'fasta')
        for i in range(aln.get_alignment_length()):
            count = Counter(aln[:, i]).most_common(1)
            # logging.info(count)
            aa_pos = i + 1
            for j in count:
                aa, freq = j
                consensus_data.append({
                    "Gene": gene,
                    "AA_cons": aa,
                    "AA_pos": aa_pos,
                    "AA_freq": freq,
                    "Sequence_type": 'Dominant',
                })
                # logging.info(f"{gene}: AA: {aa}, Position: {aa_pos}, Frequency: {freq}")

    # Convert the list of dictionaries to a DataFrame and save to CSV
    consensus_dominant_aa_df = pd.DataFrame(consensus_data)
    consensus_dominant_aa_df.to_csv(dominant_aa_path) #out_dir / 'final_core_consensus_dominant_aa_count_df.csv'
    logging.info("Finishing: preplot: find_dominant_aa")






