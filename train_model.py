import numpy as np
import os
import pandas as pd
import random
import torch
import torch.nn as nn
import torch.optim as optim

from sklearn.preprocessing import StandardScaler
from torch.utils.data import Dataset, DataLoader

conditions = ["Lupus", "Rheumatoid_Arthritis", "Multiple_Sclerosis", "Endometriosis", "Sarcoidosis", "Ulcerative_Colitis", "Crohn", "Parkinson", "Alzheimer"]
savefile = "model.pth"
train_test_split = 0.3
batch_size = 16
seed = 12345

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
    def __init__(self, num_genes, num_conditions):
        super(ConditionModel, self).__init__()
        self.input_layer = nn.Linear(num_genes, 100)
        self.layer1 = nn.Linear(100, 100)
        self.layer2 = nn.Linear(100, 100)
        self.layer3 = nn.Linear(100, num_conditions)
    
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        self.dropout = nn.Dropout(0.3)
    
    def forward(self, x):
        x = self.dropout(self.relu(self.input_layer(x)))
        x = self.dropout(self.relu(self.layer1(x)))
        x = self.dropout(self.relu(self.layer2(x)))
        x = self.sigmoid(self.layer3(x))
        return x

random.seed(seed)
np.random.seed(seed)
train_df = pd.DataFrame()
test_df = pd.DataFrame()
for condition in conditions:
    filenames = os.listdir(condition)
    random.shuffle(filenames)

    split_index = int(len(filenames) * (1 - train_test_split))
    for train_file_index in range(0, split_index):
        filename = filenames[train_file_index]
        path = f"{condition}/{filename}"
        df = pd.read_csv(path, index_col=0)
        df = np.log1p(df)
        df = df.transpose()
        df[condition] = 1
        train_df = pd.concat([train_df, df])
    
    for test_file_index in range(split_index, len(filenames)):
        filename = filenames[test_file_index]
        path = f"{condition}/{filename}"
        df = pd.read_csv(path, index_col=0)
        df = np.log1p(df)
        df = df.transpose()
        df[condition] = 1
        test_df = pd.concat([test_df, df])

train_df = train_df.fillna(0)
test_df = test_df.fillna(0)

genes = list(set(train_df.columns).difference(set(conditions)))

train_dataset = ConditionDataset(train_df, conditions)
test_dataset = ConditionDataset(test_df, conditions)

train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
test_dataloader = DataLoader(test_dataset, batch_size=batch_size, shuffle=True)

model = ConditionModel(len(train_df.columns) - len(conditions), len(conditions))
optimizer = optim.AdamW(model.parameters(), lr=0.0001, weight_decay=0.0001)
if os.path.exists(savefile):
    checkpoint = torch.load(savefile, weights_only=True)
    model.load_state_dict(checkpoint["model_state_dict"])
    optimizer.load_state_dict(checkpoint["optim_state_dict"])
model = torch.compile(model)

criterion = nn.CrossEntropyLoss()
for epoch in range(100):
    print(f"Epoch: {epoch}")
    train_metrics = {condition: {"tp": 1, "tn": 1, "fp": 1, "fn": 1} for condition in conditions}
    train_loss = 0
    for inputs, labels in train_dataloader:
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        train_loss += loss.item()
        optimizer.step()

        for batch_idx in range(len(outputs)):
            output_batch = outputs[batch_idx]
            label_batch = labels[batch_idx] 
            for condition_idx in range(len(output_batch)):
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
    train_loss /= len(train_dataset)
    print(f"Train Loss = {train_loss}")    

    print()
    model.eval()
    test_metrics = {condition: {"tp": 1, "tn": 1, "fp": 1, "fn": 1} for condition in conditions}
    test_loss = 0
    for inputs, labels in test_dataloader:
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        test_loss += loss.item()

        for batch_idx in range(len(outputs)):
            output_batch = outputs[batch_idx]
            label_batch = labels[batch_idx] 
            for condition_idx in range(len(output_batch)):
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
    test_loss /= len(test_dataset)
    print(f"Test Loss = {test_loss}")
    print()
    model.train()

    torch.save({
        "model_state_dict": model.state_dict(),
        "optim_state_dict": optimizer.state_dict(),
        "conditions": conditions
    }, savefile)