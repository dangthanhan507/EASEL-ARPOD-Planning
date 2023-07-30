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
    Q = np.eye(6)*100
    R = np.eye(3)

    pos0 = np.ones(3)*-5
    vel0 = np.ones(3)*0.1
    x0 = np.array([pos0,vel0]).reshape((6,1))
    mpc = NN_MPC(1,Q,R,horizon=5,network=network)
    xs,us = mpc.optimize(x0)
    u = us[:,0].reshape((3,1))
    
    A,B = createHCWDiscreteMatrices(1,1000,100)
    x1 = A@x0 + B@u
    x2 = A@x1 + B@us[:,1].reshape((3,1))

    out_x1 = network(torch.Tensor(np.vstack((x0,u))))
    print(xs)
    print(out_x1)
    print(x1)
    print(x2)