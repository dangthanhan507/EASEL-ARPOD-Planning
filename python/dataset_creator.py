import torch

from arpod_dynamics import createHCWMatrices, createHCWDiscreteMatrices
from controller import LQR
import numpy as np
from tqdm import tqdm
import os

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
    pos_lb, pos_hb = pos_limits
    vel_lb, vel_hb = vel_limits

    n_samples = int(n_samples)
    pos = np.random.uniform(low=pos_lb, high=pos_hb, size=(3,n_samples))
    vel = np.random.uniform(low=vel_lb, high=vel_hb, size=(3,n_samples))
    x = np.vstack((pos,vel)) #6xN
    u = np.random.uniform(low=-0.1,high=0.1,size=(3,n_samples)) #3xN
    x1 = Ak@x + Bk@u #6xN
    return np.vstack((x,u)).T, x1.T

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
    datas = torch.tensor(datas)
    labels = torch.tensor(labels)

    os.makedirs('./dataset/',exist_ok=True)
    torch.save(datas, './dataset/hcw_input.pt')
    torch.save(labels,'./dataset/hcw_labels.pt')

if __name__ == '__main__':
    load_data = True
    
    if not load_data:
        mu_GM = 1000
        R =  100
        datas, labels = createHCWDatapoints(mu_GM,R,tstep=1,n_datapoints=2e7) #allocates a few GiBs of space
        print(datas.shape, labels.shape)
        save_data(datas, labels)
    else:
        datas = torch.load('./dataset/hcw_input.pt')
        labels = torch.load('./dataset/hcw_labels.pt')
        print(datas.shape)
        print(labels.shape)
    