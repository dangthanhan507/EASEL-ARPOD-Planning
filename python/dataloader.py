import torch
from torch.utils.data import Dataset, DataLoader
import os
class DynamicsDataset(Dataset):
    def __init__(self, folder):
        
        self.datas = torch.load(os.path.join(folder, "hcw_input.pt")).to(torch.float32)
        self.labels = torch.load(os.path.join(folder, "hcw_labels.pt")).to(torch.float32)

    #overloaded operatiosn (REQUIRED for inheritance)
    def __getitem__(self, idx):
        return self.datas[idx,:], self.labels[idx,:]
    def __len__(self):
        return self.datas.shape[0]


if __name__ == '__main__':
    dataset = DynamicsDataset('./dataset/')
    dataloader = DataLoader(dataset,shuffle=True,num_workers=4,batch_size=100,pin_memory=False)
    print('Done Loading Dataset\n')