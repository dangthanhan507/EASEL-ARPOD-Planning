import torch
from arpod_dynamics import createHCWDiscreteMatrices

from controller import NN_MPC
import numpy as np
if __name__ == '__main__':

    device = 'cpu'
    network = torch.load('./models/hcw_simple_model.pt').to(device)
    Q = np.eye(6)
    R = np.eye(3)

    pos0 = np.ones(3)*5
    vel0 = np.ones(3)*0.1
    x0 = np.array([pos0,vel0]).reshape((6,1))
    mpc = NN_MPC(1,Q,R,10,network)
    xs,us = mpc.optimize(x0)
    u = us[:,0].reshape((3,1))




    tstep = 1
    mu_GM = 1000
    R =  100
    A,B = createHCWDiscreteMatrices(tstep,mu_GM,R)
    
    x1 = A@x0 + B@u

    inp = torch.Tensor(np.vstack((x0,u)))
    net_x1 = network(inp.T).T
    print(x0)
    print(x1)
    print(net_x1)