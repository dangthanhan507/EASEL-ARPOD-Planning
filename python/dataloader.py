import torch
from torch.utils.data import Dataset

from arpod_dynamics import createHCWMatrices, createHCWDiscreteMatrices
from controller import LQR
import numpy as np
from tqdm import tqdm

def runClassicARPODMission(lqr, Ak, Bk, x0):
    '''

    '''
    input_datapts = []
    label_datapts = []

    x = x0.copy()
    iter = 0
    while np.linalg.norm(x[:3,:]) > 1e-3 and iter < 1e3:
        u = lqr.getControl(x)
        x1 = Ak@x + Bk@(u)

        input_data = np.vstack((x,u))
        label_data = x1.copy()
        input_datapts.append(input_data)
        label_datapts.append(label_data)
        
        x = x1
        iter += 1

    return np.array(input_datapts).squeeze(), np.array(label_datapts).squeeze()

def runRandomHCW(Ak, Bk, pos_limits = (-10,10), vel_limits= (-1,1), n_samples=1e5):
    '''

    '''
    input_datapts = []
    label_datapts = []

    pos_lb, pos_hb = pos_limits
    vel_lb, vel_hb = vel_limits

    
    iter = 0
    while iter < n_samples:
        pos = np.random.uniform(low=pos_lb, high=pos_hb, size=(3,1))
        vel = np.random.uniform(low=vel_lb, high=vel_hb, size=(3,1))
        u = np.random.uniform(low=-0.1,high=0.1,size=(3,1))
        x = np.vstack((pos,vel))

        x1 = Ak@x + Bk@(u)

        input_data = np.vstack((x,u))
        label_data = x1.copy()
        input_datapts.append(input_data)
        label_datapts.append(label_data)
        
        iter += 1

    return np.array(input_datapts).squeeze(), np.array(label_datapts).squeeze()

def createHCWDatapoints(mu_GM, R, tstep, n_datapoints=1e6, pos_limits = (-10,10), vel_limits = (-1,1)):
    Ak,Bk = createHCWDiscreteMatrices(tstep, mu_GM, R)
    datas,labels = runRandomHCW(Ak,Bk,pos_limits,vel_limits,n_samples=n_datapoints)
    return datas, labels


def createArpodDatapoints(mu_GM, R, tstep,n_simulations=1000, pos_limits = (-10,10), vel_limits = (-1,1)):
    '''
        Description:
        ------------

        
        Parameters:
        -----------

        Returns:
        --------

    '''
    A,B = createHCWMatrices(tstep,mu_GM,R)
    Ak,Bk = createHCWDiscreteMatrices(tstep, mu_GM, R)

    Q = np.eye(6)
    R = np.eye(3)
    controller = LQR(A,B,Q,R)

    pos_lb, pos_hb = pos_limits
    vel_lb, vel_hb = vel_limits

    datas = []
    labels = []
    #run n_simulations simulations
    print("Creating Data......")
    for sim in tqdm(range(n_simulations)):
        
        pos = np.random.uniform(low=pos_lb, high=pos_hb, size=(3,1))
        vel = np.random.uniform(low=vel_lb, high=vel_hb, size=(3,1))
        x0 = np.vstack((pos,vel))
        xs,ys = runClassicARPODMission(controller, Ak, Bk, x0)

        datas.append(xs)
        labels.append(ys)
    print("Finished Creating Data")
    datas = np.vstack(datas)
    labels = np.vstack(labels)
    print(f"Data Input Shape: {datas.shape}")
    print(f"Label Shape: {labels.shape}")

    return datas, labels

def save_data(datas, labels):
    


class DynamicsDataset(Dataset):
    def __init__(self, A, B):
        pass

    #overloaded operatiosn (REQUIRED for inheritance)
    def __getitem__(self, idx):
        pass
    def __len__(self):
        pass


if __name__ == '__main__':
    mu_GM = 1000
    R =  100
    datas, labels = createHCWDatapoints(mu_GM,R,1)
    