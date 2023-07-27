import torch
from torch.utils.data import Dataset
import os
class DynamicsDataset(Dataset):
    def __init__(self, folder):
        
        self.datas = torch.load(os.path.join(folder, "hcw_input.pt"))
        self.labels = torch.load(os.path.join(folder, "hcw_labels.pt"))

    #overloaded operatiosn (REQUIRED for inheritance)
    def __getitem__(self, idx):
        return self.datas[idx,:], self.labels[idx,:]
    def __len__(self):
        return self.datas.shape[0]


if __name__ == '__main__':
    pass