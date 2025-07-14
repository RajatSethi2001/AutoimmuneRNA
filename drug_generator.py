import gymnasium as gym
import mygene
import numpy as np
import os
import pandas as pd
import random
import selfies as sf
import time
import torch
import torch.nn.functional as F
from gymnasium import spaces
from scipy import spatial
from stable_baselines3 import PPO
from train_gctx import GenePertModel, robust_minmax

def process_gene_csv(path, desired_genes, ensembl_to_gene_sym):
    df = pd.read_csv(path, index_col=0)
    df = np.log1p(df)
    df = df.apply(robust_minmax)
    df = df.transpose()

    df.columns = ensembl_to_gene_sym

    desired_genes_set = set(desired_genes)
    current_genes_set = set(df.columns)
    missing_genes = [g for g in desired_genes_set if g not in current_genes_set]
    missing_df = pd.DataFrame(
        data=0.0,
        index=df.index,
        columns=missing_genes
    )
    df = pd.concat([df, missing_df], axis=1)
    df = df.loc[:, ~df.columns.duplicated()]
    df = df[desired_genes]
    gene_expr = df.to_numpy().flatten()
    return gene_expr

def smiles_to_selfies_encoding(smiles, selfies_alphabet):
    vocab_dict = {token: i for i, token in enumerate(selfies_alphabet)}
    selfies = sf.encoder(smiles)
    tokens = list(sf.split_selfies(selfies))
    indices = [vocab_dict[token] for token in tokens]
    selfies_one_hot = np.eye(len(vocab_dict))[indices]
    return selfies_one_hot

def selfies_encoding_to_smiles(selfies_encoding, selfies_alphabet):
    selfies_argmax = np.argmax(selfies_encoding, axis=1)
    selfies = ""
    for token_idx in selfies_argmax: 
        token_value = selfies_alphabet[token_idx]
        if token_value == "[STOP]":
            break

        selfies += token_value
    
    smiles = sf.decoder(selfies)
    return smiles

class DrugGenEnv(gym.Env):
    def __init__(self, gctx_model_savefile, condition_dirs, healthy_dir, max_selfies_len=50):
        super().__init__()
        checkpoint = torch.load(gctx_model_savefile, weights_only=True)
        self.genes = checkpoint["genes"]
        self.selfies_alphabet = checkpoint["selfies_alphabet"]
        self.gctx_model = GenePertModel(len(self.genes), len(self.selfies_alphabet))
        self.gctx_model.load_state_dict(checkpoint["model_state_dict"])
        self.gctx_model.eval()
        self.max_selfies_len = max_selfies_len

        healthy_files = os.listdir(healthy_dir)
        df = pd.read_csv(f"{healthy_dir}/{healthy_files[0]}", index_col=0)
        ensembl_ids = [id.split(".")[0] for id in df.index]

        print("Converting ENSEMBL IDs to Gene Symbols")
        mg = mygene.MyGeneInfo()
        query_result = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
        id_to_symbol = {item['query']: item.get('symbol', "UNKNOWN") for item in query_result}
        ensembl_to_gene_sym = list(id_to_symbol.values())
        print("Finished Converting ENSEMBL IDs to Gene Symbols")

        self.condition_expr = {}
        for dir in condition_dirs:
            for file in os.listdir(dir):
                filename = f"{dir}/{file}"
                print(f"Processing {filename}")
                self.condition_expr[filename] = process_gene_csv(filename, self.genes, ensembl_to_gene_sym)
        
        self.healthy_expr = []
        for file in os.listdir(healthy_dir):
            filename = f"{healthy_dir}/{file}"
            print(f"Processing {filename}")
            self.healthy_expr.append(process_gene_csv(filename, self.genes, ensembl_to_gene_sym))
        
        self.observation_space = spaces.Box(low=0, high=1, shape=(len(self.genes),), dtype=np.float32)
        self.action_space = spaces.Box(low=0, high=1, shape=(self.max_selfies_len * len(self.selfies_alphabet) + 2,), dtype=np.float32)

    def step(self, action: np.ndarray):
        selfies_encoding = action[:(len(action) - 2)].reshape(self.max_selfies_len, len(self.selfies_alphabet))
        dosage_conc = action[len(action) - 2] * 2
        dosage_time = action[len(action) - 1] * 2

        with torch.no_grad():
            current_obs_expr_tensor = torch.tensor(self.current_obs_expr, dtype=torch.float32)
            selfies_encoding_tensor = torch.tensor(selfies_encoding, dtype=torch.float32)
            dosage_conc_tensor = torch.tensor(dosage_conc, dtype=torch.float32).unsqueeze(0)
            dosage_time_tensor = torch.tensor(dosage_time, dtype=torch.float32).unsqueeze(0)

            # print(current_obs_expr_tensor.shape)
            # print(selfies_encoding_tensor.shape)
            # print(dosage_conc_tensor.shape)
            # print(dosage_time_tensor.shape)
            
            selfies_encoding_tensor = F.gumbel_softmax(selfies_encoding_tensor, tau=1.0, hard=True, dim=-1)
            new_expr = self.gctx_model(current_obs_expr_tensor, selfies_encoding_tensor, dosage_conc_tensor, dosage_time_tensor).numpy().flatten()

        cos_sim_avg = 0
        for healthy_expr in self.healthy_expr:
            cos_sim_avg += (1 - spatial.distance.cosine(new_expr, healthy_expr))
        
        cos_sim_avg /= len(self.healthy_expr)
        smiles = selfies_encoding_to_smiles(selfies_encoding, self.selfies_alphabet)
        if smiles is None or smiles.strip() == "":
            cos_sim_avg = -1
        
        unnorm_dosage_conc = (10 ** dosage_conc) - 1
        unnorm_dosage_time = (10 ** dosage_time) - 1
        print(f"File: {self.current_obs_file}, SMILES = {smiles}, Dosage = {unnorm_dosage_conc} uM, Time = {unnorm_dosage_time} h, Reward = {cos_sim_avg}")

        return self.current_obs_expr, cos_sim_avg, True, False, {}

    def reset(self, seed=None, options=None):
        self.current_obs_file = random.choice(list(self.condition_expr.keys()))
        self.current_obs_expr = self.condition_expr[self.current_obs_file]
        return self.current_obs_expr, {}

gctx_model_savefile = "gctx.pth"
condition_dirs = ["Amyotrophic_Lateral_Sclerosis"]
healthy_dir = "Healthy"
policy_savefile = "test_rl"

env = DrugGenEnv(gctx_model_savefile, condition_dirs, healthy_dir)

model = PPO("MlpPolicy", env)
for epoch in range(100):
    model.learn(total_timesteps=1000, progress_bar=True)
    model.save(policy_savefile)
