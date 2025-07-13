import gymnasium as gym
import mygene
import numpy as np
import os
import pandas as pd
import torch
import torch.nn.functional as F
from gymnasium import spaces
from stable_baselines3 import PPO
from train_gctx import GenePertModel, robust_minmax

def process_gene_csv(path, desired_genes, ensembl_to_gene_sym):
    df = pd.read_csv(path, index_col=0)
    df = np.log1p(df)
    df = df.apply(robust_minmax)
    df = df.transpose()

    df.columns = ensembl_to_gene_sym
    valid_desired_genes = [g for g in desired_genes if g in df.columns]
    df = df[valid_desired_genes]

    return df.to_numpy()

class DrugGenEnv(gym.Env):
    def __init__(self, gctx_model_savefile, conditions, healthy_dir, max_selfies_len=50):
        super().__init__()
        checkpoint = torch.load(gctx_model_savefile, weights_only=True)
        self.genes = checkpoint["genes"]
        self.selfies_alphabet = checkpoint["selfies_alphabet"]
        self.gctx_model = GenePertModel(len(self.genes), len(self.selfies_alphabet))
        self.gctx_model.load_state_dict(checkpoint["model_state_dict"])
        self.max_selfies_len = max_selfies_len

        self.observation_space = spaces.Box(low=0, high=1, shape=(len(self.genes),), dtype=np.float32)
        self.action_space = spaces.Box(low=0, high=1, shape=(self.max_selfies_len * len(self.selfies_alphabet) + 2,), dtype=np.float32)

    def step(self, action: np.ndarray):
        selfies_encoding = action[:(len(action) - 2)].reshape(self.max_selfies_len, len(self.selfies_alphabet))
        selfies_encoding = F.gumbel_softmax(selfies_encoding, tau=1.0, hard=True, dim=-1)

        dosage = action[len(action) - 2] * 2
        time = action[len(action) - 1] * 2

    def reset(self):
        pass

healthy_dir = "Healthy"
healthy_files = os.listdir(healthy_dir)
df = pd.read_csv(f"{healthy_dir}/{healthy_files[0]}", index_col=0)
ensembl_ids = [id.split(".")[0] for id in df.index]

mg = mygene.MyGeneInfo()
query_result = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
id_to_symbol = {item['query']: item.get('symbol', "UNKNOWN") for item in query_result}
ensembl_to_gene_sym = list(id_to_symbol.values())

checkpoint = torch.load("gctx.pth", weights_only=True)
desired_genes = checkpoint["genes"]
print(process_gene_csv("Ulcerative_Colitis/SRR25927215.csv", desired_genes, ensembl_to_gene_sym))