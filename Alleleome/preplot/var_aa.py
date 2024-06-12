import pandas as pd

# import os
# import random
# import pathlib
# from os.path import join, getsize
# from Bio import AlignIO
# from pathlib import Path
# from Bio.Align import AlignInfo
# from Bio import SeqIO
# from Bio.Blast import NCBIWWW
# import stat
# from Bio.Blast import NCBIXML
# from Bio import SearchIO
# import numpy as np
# import timeit
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.Align.Applications import MafftCommandline
# import subprocess
# from io import StringIO
# from Bio import AlignIO
# import collections
# from collections import Counter
# import sys
# import logging
# import datetime
# import argparse


def find_variable_aa(aa_vars_path, variable_aa_path):
    print('Start running var_aa.py')

    df = pd.read_csv(aa_vars_path, na_filter=False)

    df.drop(index=df.loc[df.AA_mutation_type != 'Substitution'].index, inplace=True)
    df.reset_index(inplace=True)

    df['Cons_aa'] = df['Cons_aa'].apply(list)
    df['Var_aa'] = df['Var_aa'].apply(list)
    s = pd.Series([np.arange(x.AA_start_pos, x.AA_end_pos+1, 1) for x in df.itertuples()], index=df.index)
    df['AA_pos'] = pd.Series(data=s, index=df.index)
    df_f = df.explode(['AA_pos','Cons_aa','Var_aa'])

    df_f.to_csv(variable_aa_path)


# def find_dominant_var_all(variable_aa_path, ):
#     print('Start running dominant_var_all.py')

#     #This is a file containing all separeted substitution positions along with variant details
#     df = pd.read_csv(variable_aa_path, low_memory=False)

#     df_var_1 = df[['AA_pos','Var_aa','Sequence_type']].groupby(df1['Gene']).value_counts().reset_index(name='AA_freq')

#     df_var_1.to_csv(alleleome_dir_path + 'final_pan_aa_thresh_core_var_genome_count.csv')

#     df_dom = pd.read_csv(alleleome_dir_path + 'final_core_consensus_dominant_aa_count_df.csv', low_memory=False)

#     col_dict_dom = {'AA_cons': 'Amino_acid', 'AA_freq':'Genome_count'}   ## key→old name, value→new name

#     df_dom.columns = [col_dict_dom.get(x, x) for x in df_dom.columns]
#     df_dom

#     col_dict_var = {'Var_aa': 'Amino_acid', 'AA_freq':'Genome_count'}   ## key→old name, value→new name

#     df_var_1.columns = [col_dict_var.get(x, x) for x in df_var_1.columns]

#     df_var_1
#     df_dom_var = pd.concat([df_var_1, df_dom]).sort_values("Gene")

#     df_dom_var.to_csv(alleleome_dir_path + 'final_pan_aa_thresh_core_genes_dominant_variant_genome_count_pos.csv')

#     df_dom_var=pd.read_csv(alleleome_dir_path + 'final_pan_aa_thresh_core_genes_dominant_variant_genome_count_pos.csv', low_memory=False)

#     df_gap=df_dom_var.loc[(df_dom_var["Amino_acid"] == '-') , ["Gene","AA_pos"]]


#     df_gap.to_csv(alleleome_dir_path + 'pan_aa_thresh_core_genes_aa_pos_with_gaps.csv')

#     df_dom_var=pd.merge(df_dom_var, df_gap, on=["Gene", "AA_pos"], how="outer", indicator=True)

#     df_dom_var.to_csv(alleleome_dir_path + 'pan_aa_thresh_core_genes_aa_pos_with_gaps_in_dom_and_var.csv')
#     df_filt = df_dom_var.loc[df_dom_var["_merge"] == "left_only"].drop("_merge", axis=1)
#     df_filt.to_csv(alleleome_dir_path + 'pan_aa_thresh_core_genes_dom_var_filtered_pos_with_gaps.csv')

#     df2=pd.read_csv(alleleome_dir_path + 'pan_aa_thresh_core_genes_dom_var_filtered_pos_with_gaps.csv', index_col=0, low_memory=False)

#     cols = ['Genome_count']

#     g = df2.groupby('Gene')[cols]
#     min1 = g.transform('min')

#     max1 = g.transform('max')

#     df1 = df2.join(df2[cols].sub(min1).div(max1 - min1).add_suffix('_norm'))

#     df1['Genome_count_norm']=df1['Genome_count_norm'].apply(lambda x:round(x,2))
#     df1.to_csv(alleleome_dir_path + 'final_pan_aa_thresh_core_genes_dom_var_genome_count_pos_normalized.csv')
    
#     df_summary = pd.read_csv(os.path.join(alleleome_dir_path, 'df_pangene_summary_v2.csv'), low_memory=False)
#     core_gene_list = list(df_summary[df_summary['pangenome_class_2'] == 'Core']['Gene'])

#     #for getting the single gene details
#     for core_gene_name in core_gene_list: 
#         if not os.path.exists(os.path.join(output_folder, species, 'alleleome', core_gene_name)):
#             os.makedirs(os.path.join(output_folder, species, 'alleleome', core_gene_name))

#         df2[df2.Gene == core_gene_name].to_csv(os.path.join(output_folder, species,'alleleome', core_gene_name, core_gene_name + '_pan_aa_thresh_core_dom_var_pos.csv'))






