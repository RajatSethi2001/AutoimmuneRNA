import numpy as np
import os
import pandas as pd
import random
from sklearn.feature_selection import f_classif
from sklearn.preprocessing import LabelEncoder
from utils import get_zscores, get_zscore_minmax

def main():
    gene_count = 5000
    condition_parent_dir = "Conditions"
    directory_list = ["Lupus", "Podoconiosis", "Primary_Sclerosing_Cholangitis", "Ulcerative_Colitis", "Shingles", "Sepsis", "Scleroderma", "MRSA_Bacteremia", "Crohns_Disease", "Acute_Pancreatitis", "Aneurysm", "Tuberculosis", "Acute_Myeloid_Leukemia", "Endocarditis", "Schistosomiasis", "Leprosy", "Amyotrophic_Lateral_Sclerosis", "Chronic_Myeloid_Leukemia", "Dengue", "Alzheimer", "Coronary_Artery_Disease", "COPD", "Breast_Cancer", "Crimean_Congo_Hemorrhagic_Fever"]
    
    gene_expr_df = pd.DataFrame()
    labels = []
    for directory in directory_list:
        path = f"{condition_parent_dir}/{directory}"
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
            gene_expr_df = pd.concat([gene_expr_df, df])
            labels.append(directory)

    label_encoder = LabelEncoder()
    labels_encoded = label_encoder.fit_transform(labels)
    f, pvals = f_classif(gene_expr_df.values, labels_encoded)
    gene_ranking = pd.Series(f, index=gene_expr_df.columns).sort_values(ascending=False)
    top_genes = gene_ranking.head(gene_count).index.tolist()

    with open("Data/important_genes.txt", "w") as f:
        for gene_id in top_genes:
            f.write(f"{gene_id}\n")

if __name__=="__main__":
    main()