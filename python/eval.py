import torch
from arpod_dynamics import createHCWDiscreteMatrices
from networks import SimpleDynamicsNetwork
from controller import NN_MPC
import numpy as np
if __name__ == '__main__':

    device = 'cpu'
    network = SimpleDynamicsNetwork()
    network_load = torch.load('./models/hcw_simple_model.pt').to(device)
    network.load_state_dict(network_load.state_dict())
    A,B = createHCWDiscreteMatrices(1,1000,100)

    Q = np.eye(6)*100
    R = np.eye(3)

    pos0 = np.ones(3)*-5
    vel0 = np.ones(3)*0.1
    x0 = np.array([pos0,vel0]).reshape((6,1))
    mpc = NN_MPC(1,Q,R,horizon=5,network=network)
    
    t = 0
    finalT = 100
    while t < finalT:
        xs,us = mpc.optimize(x0)
        u = us[:,0].reshape((3,1))
        x0 = A@x0 + B@u

        t += 1

    print(x0)
    

    