import pandas as pd
import numpy as np
from pathlib import Path


def find_variable_aa(aa_vars_path, variable_aa_path):
    print("Start running var_aa.py")

    df = pd.read_csv(aa_vars_path, na_filter=False)

    df.drop(index=df.loc[df.AA_mutation_type != "Substitution"].index, inplace=True)
    df.reset_index(inplace=True)

    df["Cons_aa"] = df["AA_cons_seq"].apply(list)
    df["Var_aa"] = df["AA_seq_change"].apply(list)
    s = pd.Series(
        [np.arange(x.AA_start_pos, x.AA_end_pos + 1, 1) for x in df.itertuples()],
        index=df.index,
    )
    df["AA_pos"] = pd.Series(data=s, index=df.index)
    df_f = df.explode(["AA_pos", "Cons_aa", "Var_aa"])

    df_f.to_csv(variable_aa_path)


def find_dominant_var_all(
    variable_aa_path,
    dominant_aa_path,
    dom_var_path,
    gaps_path,
    filt_norm_path,
    dom_var_out_dir,
):
    print("Start running dominant_var_all.py")
    dom_var_out_dir = Path(dom_var_out_dir)

    # This is a file containing all separeted substitution positions along with variant details
    df = pd.read_csv(variable_aa_path, low_memory=False)

    df_var = (
        df[["AA_pos", "Var_aa", "Sequence_type"]]
        .groupby(df["Gene"])
        .value_counts()
        .reset_index(name="AA_freq")
    )
    # df_dom = pd.read_csv(
    #     alleleome_dir_path + "final_core_consensus_dominant_aa_count_df.csv",
    #     low_memory=False,
    # )

    df_dom.rename(
        columns={"AA_cons": "Amino_acid", "AA_freq": "Genome_count"}, inplace=True
    )
    df_var.rename(
        columns={"Var_aa": "Amino_acid", "AA_freq": "Genome_count"}, inplace=True
    )

    df_dom_var = pd.concat([df_var, df_dom]).sort_values("Gene")

    df_dom_var.to_csv(dom_var_path)

    df_gap = df_dom_var.loc[(df_dom_var["Amino_acid"] == "-"), ["Gene", "AA_pos"]]
    df_gap.to_csv(gaps_path)

    df_filt = df_dom_var.loc[(df_dom_var["Amino_acid"] != "-"), :]
    cols = ["Genome_count"]
    g = df_filt.groupby("Gene")[cols]
    min1 = g.transform("min")
    max1 = g.transform("max")

    df_norm = df_filt.join(df_filt[cols].sub(min1).div(max1 - min1).add_suffix("_norm"))

    df_norm["Genome_count_norm"] = df_norm["Genome_count_norm"].apply(
        lambda x: round(x, 2)
    )
    df_norm.to_csv(filt_norm_path)

    # for getting the single gene details
    for gene in gene_list:
        gene_path = dom_var_out_dir / gene
        gene_path.mkdir(parents=True, exist_ok=True)

        df_norm[df_norm.Gene == gene].to_csv(
            gene_path / f"{gene}_pan_aa_thresh_core_dom_var_pos.csv"
        )
