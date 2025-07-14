import h5py
import numpy as np
import os
import pandas as pd
import random
import selfies as sf
import torch
import torch.nn as nn
import torch.optim as optim
from cmapPy.pandasGEXpress import parse
from torch.utils.data import Dataset, DataLoader, Subset
"""
/0
/0/DATA
/0/DATA/0
/0/DATA/0/matrix
/0/META
/0/META/COL
/0/META/COL/cell_id
/0/META/COL/distil_id
/0/META/COL/id
/0/META/COL/pert_dose
/0/META/COL/pert_dose_unit
/0/META/COL/pert_id
/0/META/COL/pert_idose
/0/META/COL/pert_iname
/0/META/COL/pert_itime
/0/META/COL/pert_time
/0/META/COL/pert_time_unit
/0/META/COL/pert_type
/0/META/ROW
/0/META/ROW/id
/0/META/ROW/pr_gene_symbol
/0/META/ROW/pr_gene_title
/0/META/ROW/pr_is_bing
/0/META/ROW/pr_is_lm
"""
def clean_dose_unit(raw):
    if isinstance(raw, bytes):
        if raw in [b'\xfd\xfdM', b'\xc2\xb5M']:
            return 'uM'
        elif raw == b'\xfd\xfdL':
            return 'uL'
        elif raw == b'ng/mL':
            return 'ng/mL'
        elif raw == b'ng':
            return 'ng'
        elif raw == b'ng/\xfd\xfdL':
            return 'ng/uL'
        elif raw == b'-666':
            return None
        else:
            return raw.decode('utf-8', errors='ignore')  # fallback
    return raw

def robust_minmax(vector, lower=5, upper=95):
    lo = np.percentile(vector, lower)
    hi = np.percentile(vector, upper)
    return np.clip((vector - lo) / (hi - lo), 0, 1)

def encode_smiles(smiles, selfies_alphabet):
    vocab_dict = {token: i for i, token in enumerate(selfies_alphabet)}
    selfies = sf.encoder(smiles)
    tokens = list(sf.split_selfies(selfies))
    indices = [vocab_dict[token] for token in tokens]
    selfies_one_hot = np.eye(len(vocab_dict))[indices]
    return selfies_one_hot

class GenePertDataset(Dataset):
    def __init__(self, gctx_file, compound_file):
        self.gctx_fp = h5py.File(gctx_file, "r")
        self.cell_types = [s.decode('utf-8') for s in self.gctx_fp["0/META/COL/cell_id"]]
        self.pert_dose = [float(s.decode('utf-8').split("|")[0]) for s in self.gctx_fp["0/META/COL/pert_dose"]]
        self.pert_dose_units = [clean_dose_unit(s) for s in self.gctx_fp["0/META/COL/pert_dose_unit"]]
        self.pert_time = [float(s) for s in self.gctx_fp["0/META/COL/pert_time"]]
        self.pert_time_units = [s.decode('utf-8') for s in self.gctx_fp["0/META/COL/pert_time_unit"]]
        self.pert_types = [s.decode('utf-8') for s in self.gctx_fp["0/META/COL/pert_type"]]
        self.pert_id = [s.decode('utf-8') for s in self.gctx_fp["0/META/COL/pert_id"]]
        self.gene_symbols = [s.decode("utf-8") for s in self.gctx_fp["/0/META/ROW/pr_gene_symbol"]]
        
        self.ctl_idx = [i for i, val in enumerate(self.pert_types) if val == "ctl_untrt" or val == "ctl_vehicle"]
        random.shuffle(self.ctl_idx)
        self.ctl_idx = self.ctl_idx[:2000]
        self.trt_idx = [i for i, val in enumerate(self.pert_types) if val == "trt_cp"]

        self.valid_trt_samples = {}
        for idx in self.trt_idx:
            if self.pert_dose_units[idx] == "uM" and self.pert_time_units[idx] == "h":
                cell_type = self.cell_types[idx]
                if cell_type in self.valid_trt_samples:
                    self.valid_trt_samples[cell_type].append(idx)
                else:
                    self.valid_trt_samples[cell_type] = [idx]
        
        compound_df = pd.read_csv(compound_file, sep="\t")
        self.smiles_lookup = compound_df.set_index("pert_id")["canonical_smiles"].to_dict()
        selfies_list = []
        for cell_type, trt_indices in self.valid_trt_samples.items():
            idx = 0
            while idx < len(trt_indices):
                trt_idx = trt_indices[idx]
                try:
                    smiles = self.smiles_lookup[self.pert_id[trt_idx]]
                    selfies = sf.encoder(smiles)
                    selfies_list.append(selfies)
                except:
                    trt_indices.remove(trt_idx)
                else:
                    idx += 1

        self.selfies_alphabet = list(sf.get_alphabet_from_selfies(selfies_list))
        self.selfies_alphabet.append(".")
        self.selfies_alphabet.append("[STOP]")
        self.num_genes = self.gctx_fp["0/DATA/0/matrix"].shape[1]

        idx = 0
        while idx < len(self.ctl_idx):
            ctl_idx = self.ctl_idx[idx]
            cell_type = self.cell_types[ctl_idx]
            if cell_type in self.valid_trt_samples and len(self.valid_trt_samples[cell_type]) > 0:
                idx += 1
            else:
                self.ctl_idx.pop(idx)

        del self.pert_dose_units
        del self.pert_time_units
        del self.pert_types

    def __len__(self):
        return len(self.ctl_idx)

    def __getitem__(self, idx):
        ctl_idx = self.ctl_idx[idx]
        cell_type = self.cell_types[ctl_idx]

        trt_idx = random.choice(self.valid_trt_samples[cell_type])
        dose = self.pert_dose[trt_idx]
        time = self.pert_time[trt_idx]
        smiles = self.smiles_lookup[self.pert_id[trt_idx]]
        selfies_one_hot = encode_smiles(smiles, self.selfies_alphabet)

        ctl_expr = self.gctx_fp["0/DATA/0/matrix"][ctl_idx, :]
        trt_expr = self.gctx_fp["0/DATA/0/matrix"][trt_idx, :]
        
        ctl_expr = torch.tensor(robust_minmax(ctl_expr), dtype=torch.float32)
        trt_expr = torch.tensor(robust_minmax(trt_expr), dtype=torch.float32)
        selfies_one_hot = torch.tensor(selfies_one_hot, dtype=torch.float32)
        dose = torch.tensor(np.log10(dose + 1), dtype=torch.float32).unsqueeze(0)
        time = torch.tensor(np.log10(time + 1), dtype=torch.float32).unsqueeze(0)

        return ctl_expr, trt_expr, selfies_one_hot, dose, time

    def get_selfies_alphabet(self):
        return self.selfies_alphabet

    def get_gene_symbols(self):
        return self.gene_symbols

class GenePertModel(nn.Module):
    def __init__(self, num_genes, alphabet_len, ctl_fc_size=200):
        super().__init__()
        self.selfies_rnn = nn.GRU(
            input_size=alphabet_len,
            hidden_size=alphabet_len,
            num_layers=3,
            batch_first=True,
            dropout=0.2
        )

        self.ctl_fc_size = ctl_fc_size
        self.ctl_fc = nn.Linear(num_genes, ctl_fc_size)

        self.fc_input_len = ctl_fc_size + alphabet_len + 2
        self.fc = nn.Linear(self.fc_input_len, self.fc_input_len)
        self.output = nn.Linear(self.fc_input_len, num_genes)

        self.activation = nn.GELU()
    
    def forward(self, ctl_expr, selfies_one_hot, dose, time):
        _, h_n = self.selfies_rnn(selfies_one_hot)
        h_n = self.activation(h_n[-1])
        ctl_fc = self.activation(self.ctl_fc(ctl_expr))
        fc_input = torch.cat((ctl_fc, h_n, dose, time), dim=-1)
        fc = self.activation(self.fc(fc_input))
        output = self.activation(self.output(fc))
        return output

def main():
    seed = 100
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    train_test_split = 0.2
    dataset = GenePertDataset("annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx", "compoundinfo_beta.txt")
    model_savefile = "gctx.pth"

    train_size = int(len(dataset) * (1 - train_test_split))
    test_size = int(len(dataset) * train_test_split)
    indices = list(range(len(dataset)))
    random.shuffle(indices)
    train_indices = indices[:train_size]
    test_indices = indices[train_size:]

    train_dataset = Subset(dataset, train_indices)
    test_dataset = Subset(dataset, test_indices)

    batch_size = 1
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=True)

    criterion = nn.MSELoss()
    model = GenePertModel(len(dataset.get_gene_symbols()), len(dataset.get_selfies_alphabet()))
    optimizer = optim.AdamW(model.parameters(), lr=0.001, weight_decay=0.0001)

    if os.path.exists(model_savefile):
        checkpoint = torch.load(model_savefile, weights_only=True)
        model.load_state_dict(checkpoint["model_state_dict"])
        optimizer.load_state_dict(checkpoint["optim_state_dict"])

    epochs = 100

    for epoch in range(epochs):
        print(f"Training Epoch {epoch}")
        model.train()
        train_loss = 0
        batch = 0
        for ctl_expr, trt_expr, selfies_one_hot, dose, time in train_loader:
            optimizer.zero_grad()
            output = model(ctl_expr, selfies_one_hot, dose, time)
            loss = criterion(output, trt_expr)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()
            batch += 1

            if batch % 100 == 0:
                print(f"Train Batch {batch}: Loss = {loss.item()}")
        
        train_loss /= train_size
        print(f"Training Loss = {train_loss}")
        
        print(f"Saving Model Data: {model_savefile}")
        torch.save({
            "model_state_dict": model.state_dict(),
            "optim_state_dict": optimizer.state_dict(),
            "genes": dataset.get_gene_symbols(),
            "selfies_alphabet": dataset.get_selfies_alphabet()
        }, model_savefile)

        print(f"Testing Epoch {epoch}")
        model.eval()
        test_loss = 0
        batch = 0
        for ctl_expr, trt_expr, selfies_one_hot, dose, time in test_loader:
            output = model(ctl_expr, selfies_one_hot, dose, time)
            loss = criterion(output, trt_expr)
            test_loss += loss.item()

            batch += 1

            if batch % 100 == 0:
                print(f"Test Batch {batch}: Loss = {loss.item()}")

        test_loss /= test_size    
        print(f"Testing Batch Loss = {test_loss}")

if __name__=="__main__":
    main()