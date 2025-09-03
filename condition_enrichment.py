import gseapy as gp
import mygene
import numpy as np
import os
import pandas as pd
import random
import torch
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from utils import get_zscores, get_zscore_minmax

def main():
    with open("Data/important_genes.txt", "r") as f:
        gene_ids = f.read().splitlines()

    condition_parent_dir = "Conditions"
    condition_of_interest = "Hidradenitis_Supparativa"
    other_conditions = ["Amyotrophic_Lateral_Sclerosis", "Lupus", "Parkinson", "Podoconiosis", "Primary_Sclerosing_Cholangitis", "Ulcerative_Colitis", "Shingles", "Sepsis", "Scleroderma", "MRSA_Bacteremia", "Crohns_Disease", "Acute_Pancreatitis", "Aneurysm", "Tuberculosis", "Acute_Myeloid_Leukemia", "Endocarditis", "Schistosomiasis", "Leprosy", "Chronic_Myeloid_Leukemia", "Dengue", "Alzheimer", "Coronary_Artery_Disease", "COPD", "Breast_Cancer", "Crimean_Congo_Hemorrhagic_Fever", "Hypertension,Drug_Abuse", "Hypertension", "COVID19", "Depression", "PTSD", "HIV", "HIV,Tuberculosis", "Malaria", "SFTS", "Cystic_Fibrosis", "Chikungunya", "Rheumatoid_Arthritis", "Polycystic_Kidney_Disease", "Myelofibrosis"]

    condition_df = pd.DataFrame()
    path = f"{condition_parent_dir}/{condition_of_interest}"
    print(f"Processing Directory: {path}")
    filenames = os.listdir(path)
    random.shuffle(filenames)

    for filename in filenames:
        csv_path = f"{path}/{filename}"
        df = pd.read_csv(csv_path, index_col=0)
        df = np.log2(df + 1)
        df = df.apply(get_zscores, axis=0)
        df = df.apply(get_zscore_minmax, axis=0)
        df = df.transpose()
        df = df.loc[:, gene_ids]
        condition_df = pd.concat([condition_df, df])

    other_condition_df = pd.DataFrame()
    for other_condition in other_conditions:
        path = f"{condition_parent_dir}/{other_condition}"
        print(f"Processing Directory: {path}")
        filenames = os.listdir(path)
        random.shuffle(filenames)

        for filename in filenames:
            csv_path = f"{path}/{filename}"
            df = pd.read_csv(csv_path, index_col=0)
            df = np.log2(df + 1)
            df = df.apply(get_zscores, axis=0)
            df = df.apply(get_zscore_minmax, axis=0)
            df = df.transpose()
            df = df.loc[:, gene_ids]
            other_condition_df = pd.concat([other_condition_df, df])

    gene_pvals = []
    for gene in gene_ids:
        condition_vals = condition_df[gene].values
        other_condition_vals = other_condition_df[gene].values
        stat, pval = mannwhitneyu(condition_vals, other_condition_vals, alternative="two-sided")
        gene_pvals.append(pval)

    gene_adj_pvals = multipletests(gene_pvals, method="fdr_bh")[1]
    mean_diff = condition_df.mean() - other_condition_df.mean()
    gene_scores = -np.log10(gene_adj_pvals) * np.sign(mean_diff)

    ensembl_ids = [x.split(".")[0] for x in gene_ids]
    mg = mygene.MyGeneInfo()
    query_result = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
    id_to_symbol = {item['query']: item.get('symbol', None) for item in query_result}
    gene_symbols = list(filter(None, id_to_symbol.values()))

    ranked_genes = dict(zip(gene_symbols, gene_scores))
    rnk = pd.Series(ranked_genes).sort_values(ascending=False) 

    gsea_res = gp.prerank(rnk=rnk,
                        gene_sets=['KEGG_2021_Human', 'Reactome_2022', 'MSigDB_Hallmark_2020', 'GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021'],
                        outdir=f'Enrichments/{condition_of_interest}')

    print(gsea_res.res2d.head())

if __name__=="__main__":
    main()