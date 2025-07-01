import numpy as np
import os
import pandas as pd
import random
import time
import torch
import torch.nn as nn
import torch.optim as optim

from sklearn.preprocessing import StandardScaler
from torch.utils.data import Dataset, DataLoader

directory_list = ["Lupus", "Podoconiosis", "Primary_Sclerosing_Cholangitis", "Ulcerative_Colitis", "Shingles", "Sepsis", "Scleroderma", "MRSA_Bacteremia", "Crohns_Disease", "Acute_Pancreatitis", "Aneurysm", "Tuberculosis", "Acute_Myeloid_Leukemia", "Endocarditis", "Schistosomiasis", "Leprosy", "Amyotrophic_Lateral_Sclerosis", "Chronic_Myeloid_Leukemia", "Dengue", "Alzheimer", "Restless_Legs_Syndrome", "Coronary_Artery_Disease", "COPD", "Breast_Cancer", "Crimean_Congo_Hemorrhagic_Fever", "Hypertension,Drug_Abuse", "Hypertension", "COVID19", "Depression", "PTSD", "HIV", "HIV,Tuberculosis", "Malaria", "Hidradenitis_Supparativa", "SFTS", "Cystic_Fibrosis", "Chikungunya", "Rheumatoid_Arthritis"]
savefile = "model.pth"
train_test_split = 0.2
batch_size = 16
seed = 123456

class ConditionDataset(Dataset):
    def __init__(self, df: pd.DataFrame, conditions):
        self.df = df
        self.conditions = conditions
        self.genes = list(set(df.columns).difference(set(conditions)))
    
    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        gene_data = torch.tensor(self.df.iloc[idx, :][self.genes].values, dtype=torch.float32)
        condition_data = torch.tensor(self.df.iloc[idx, :][self.conditions].values, dtype=torch.float32)
        return gene_data, condition_data

class ConditionModel(nn.Module):
    def __init__(self, num_genes, num_conditions, hidden_neurons=100):
        super().__init__()
        self.input_layer = nn.Linear(num_genes, hidden_neurons)
        self.bn0 = nn.BatchNorm1d(hidden_neurons)
        self.layer1 = nn.Linear(hidden_neurons, hidden_neurons)
        self.bn1 = nn.BatchNorm1d(hidden_neurons)
        self.layer2 = nn.Linear(hidden_neurons, hidden_neurons)
        self.bn2 = nn.BatchNorm1d(hidden_neurons)
        self.output_layer = nn.Linear(hidden_neurons, num_conditions)
    
        self.activation = nn.GELU()
        self.dropout = nn.Dropout(0.4)
    
    def forward(self, x):
        x_min = x.min(dim=0, keepdim=True).values
        x_max = x.max(dim=0, keepdim=True).values
        denominator = x_max - x_min
        denominator[denominator == 0] = 1e-6
        x = (x - x_min) / denominator
        x = self.dropout(self.activation(self.bn0(self.input_layer(x))))
        x = self.dropout(self.activation(self.bn1(self.layer1(x))))
        x = self.dropout(self.activation(self.bn2(self.layer2(x))))
        x = self.output_layer(x)
        return x

random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
train_df = pd.DataFrame()
test_df = pd.DataFrame()
conditions = set()
for directory in directory_list:
    filenames = os.listdir(directory)
    random.shuffle(filenames)

    split_index = int(len(filenames) * (1 - train_test_split))
    for train_file_index in range(0, split_index):
        filename = filenames[train_file_index]
        path = f"{directory}/{filename}"
        df = pd.read_csv(path, index_col=0)
        df = np.log1p(df)
        df = df.transpose()
        dir_conditions = directory.split(",")
        for condition in dir_conditions:
            df[condition] = 1
            conditions.add(condition)
        train_df = pd.concat([train_df, df])
    
    for test_file_index in range(split_index, len(filenames)):
        filename = filenames[test_file_index]
        path = f"{directory}/{filename}"
        df = pd.read_csv(path, index_col=0)
        df = np.log1p(df)
        df = df.transpose()
        dir_conditions = directory.split(",")
        for condition in dir_conditions:
            df[condition] = 1
            conditions.add(condition)
        test_df = pd.concat([test_df, df])

train_df = train_df.fillna(0)
test_df = test_df.fillna(0)

genes = list(set(train_df.columns).difference(conditions))
conditions = list(conditions)

threshold = 0.8
zero_fraction = (train_df[genes] <= 0.5).sum(axis=0) / len(train_df)
genes = list(train_df[genes].loc[:, zero_fraction < threshold].columns)

train_df = train_df[genes + conditions]
test_df = test_df[genes + conditions]

train_dataset = ConditionDataset(train_df, conditions)
test_dataset = ConditionDataset(test_df, conditions)

train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=True)
test_dataloader = DataLoader(test_dataset, batch_size=batch_size, shuffle=True, drop_last=True)

model = ConditionModel(len(genes), len(conditions))
optimizer = optim.AdamW(model.parameters(), lr=0.001, weight_decay=0.0001)
model = torch.compile(model)
if os.path.exists(savefile):
    checkpoint = torch.load(savefile, weights_only=True)
    model.load_state_dict(checkpoint["model_state_dict"])
    optimizer.load_state_dict(checkpoint["optim_state_dict"])

criterion = nn.BCEWithLogitsLoss()
random.seed(int(time.time()))
np.random.seed(int(time.time()))
torch.manual_seed(int(time.time()))
for epoch in range(100):
    print(f"Epoch: {epoch}")
    train_metrics = {condition: {"tp": np.float32(0), "tn": np.float32(0), "fp": np.float32(0), "fn": np.float32(0)} for condition in conditions}
    train_loss = 0
    model.train()
    for inputs, labels in train_dataloader:
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        train_loss += loss.item()
        optimizer.step()

        for batch_idx in range(len(outputs)):
            output_batch = torch.sigmoid(outputs[batch_idx])
            label_batch = torch.sigmoid(labels[batch_idx]) 
            for condition_idx in range(len(conditions)):
                condition = conditions[condition_idx]
                output_choice = round(output_batch[condition_idx].item())
                label_choice = round(label_batch[condition_idx].item())

                if output_choice == 1 and label_choice == 1:
                    train_metrics[condition]["tp"] += 1
                
                elif output_choice == 0 and label_choice == 0:
                    train_metrics[condition]["tn"] += 1
                
                elif output_choice == 1 and label_choice == 0:
                    train_metrics[condition]["fp"] += 1
                
                elif output_choice == 0 and label_choice == 1:
                    train_metrics[condition]["fn"] += 1
    
    for condition in conditions:
        tp = train_metrics[condition]["tp"]
        tn = train_metrics[condition]["tn"]
        fp = train_metrics[condition]["fp"]
        fn = train_metrics[condition]["fn"]
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f_score = 2 * precision * recall / (precision + recall)
        print(f"{condition} Train Metrics: Precision = {round(precision, 3)}, Recall = {round(recall, 3)}, F-Score = {round(f_score, 3)}")
    
    train_loss /= (len(train_dataset) / batch_size)
    print(f"Train Loss = {train_loss}")    

    print()
    model.eval()
    test_metrics = {condition: {"tp": np.float32(0), "tn": np.float32(0), "fp": np.float32(0), "fn": np.float32(0)} for condition in conditions}
    test_loss = 0
    for inputs, labels in test_dataloader:
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        test_loss += loss.item()

        for batch_idx in range(len(outputs)):
            output_batch = torch.sigmoid(outputs[batch_idx])
            label_batch = torch.sigmoid(labels[batch_idx]) 
            for condition_idx in range(len(conditions)):
                condition = conditions[condition_idx]
                output_choice = round(output_batch[condition_idx].item())
                label_choice = round(label_batch[condition_idx].item())

                if output_choice == 1 and label_choice == 1:
                    test_metrics[condition]["tp"] += 1
                
                elif output_choice == 0 and label_choice == 0:
                    test_metrics[condition]["tn"] += 1
                
                elif output_choice == 1 and label_choice == 0:
                    test_metrics[condition]["fp"] += 1
                
                elif output_choice == 0 and label_choice == 1:
                    test_metrics[condition]["fn"] += 1
    
    for condition in conditions:
        tp = test_metrics[condition]["tp"]
        tn = test_metrics[condition]["tn"]
        fp = test_metrics[condition]["fp"]
        fn = test_metrics[condition]["fn"]
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f_score = 2 * precision * recall / (precision + recall)
        print(f"{condition} Test Metrics: Precision = {round(precision, 3)}, Recall = {round(recall, 3)}, F-Score = {round(f_score, 3)}")
    
    test_loss /= (len(test_dataset) / batch_size)
    print(f"Test Loss = {test_loss}")
    print()

    torch.save({
        "model_state_dict": model.state_dict(),
        "optim_state_dict": optimizer.state_dict(),
        "conditions": conditions
    }, savefile)