import pandas as pd
from pathlib import Path

def calculate_dn_ds(codon_muts_path, dn_ds_path):
    print('Start running dn_ds_count.py')
    
    df = pd.read_csv(codon_muts_path)
    df = df.groupby(['Gene','Cons_codon','Query_codon','Codon_position','AA_mutation_effect'])['GCF_id'].count().reset_index()

    # Remove rows where 'Query_codon' contains 'N'
    df_without_n = df[~df['Query_codon'].str.contains('N')]
    df_per_gene = df_without_n.groupby(['Gene','AA_mutation_effect'])['GCF_id'].sum().reset_index(name='Total_per_gene')

    df_per_gene = df_per_gene.pivot(index=['Gene'], columns='AA_mutation_effect') #.fillna(0)
    df_per_gene.columns = df_per_gene.columns.to_flat_index().map(lambda x: f"{x[1].replace('-', '_')}_count")
    df_per_gene = df_per_gene.reset_index().rename_axis(None, axis=1)

    df_per_gene['dN/dS'] = (df_per_gene['Non_synonymous_count']/df_per_gene['Synonymous_count'])

    df_per_gene = df_per_gene.dropna(axis=0)
    df_per_gene.to_csv(dn_ds_path)