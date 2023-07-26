import torch
from torch.utils.data import Dataset

from arpod_dynamics import createHCWMatrices, createHCWDiscreteMatrices
from controller import LQR


def createArpodDatapoints(data_size, box_limits):
    '''
        Description:
        ------------

        
        Parameters:
        -----------

        Returns:
        --------


    '''


    pass


class DynamicsDataset(Dataset):
    def __init__(self, A, B):
        pass

    #overloaded operatiosn (REQUIRED for inheritance)
    def __getitem__(self, idx):
        pass
    def __len__(self):
        pass