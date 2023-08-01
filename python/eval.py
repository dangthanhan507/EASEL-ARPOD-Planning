import torch
from torch.autograd import Variable
from mpc import mpc
from mpc.mpc import QuadCost, LinDx

import torch.nn as nn
class SimpleDynamicsNetwork(nn.Module):
    def __init__(self):
        super().__init__() #must call

        #setting up hidden layer

        # control input is 3, state input is 6... that makes 9
        self.layers = []
        self.layers.append(nn.Linear(9,64,bias=True))
        self.layers.append(nn.Linear(64,128,bias=True))
        self.layers.append(nn.Linear(128,256,bias=True))
        self.layers.append(nn.Linear(256,256,bias=True))
        self.layers.append(nn.Linear(256,256,bias=True))
        self.layers.append(nn.Linear(256,128,bias=True))
        self.layers.append(nn.Linear(128,64,bias=True))
        self.layers.append(nn.Linear(64,6,bias=True))
        self.model = nn.Sequential(*self.layers)
    
    def forward(self, x,u):
        '''
        '''
        inp = torch.hstack((x,u))
        return self.model(inp)

if __name__ == '__main__':
    torch.manual_seed(0)

    device = 'cpu'
    network = SimpleDynamicsNetwork()
    network_load = torch.load('./models/hcw_simple_model.pt').to(device)
    network.load_state_dict(network_load.state_dict())

    n_batch, n_state, n_ctrl, T = 1, 6, 3, 10
    n_sc = n_state + n_ctrl

    #load network

    # Randomly initialize a PSD quadratic cost
    C = torch.randn(T*n_batch, n_sc, n_sc)
    C = torch.bmm(C, C.transpose(1, 2)).view(T, n_batch, n_sc, n_sc)
    c = torch.randn(T, n_batch, n_sc)

    # The initial state.
    n_batch = 1
    x_init = torch.Tensor([[10,10,10,0,0,0]]) #1,6
    ctrl_mpc = mpc.MPC(n_state,
                        n_ctrl,
                        T,
                        u_lower=-0.1*torch.ones(T,n_batch,n_ctrl),
                        u_upper=0.1*torch.ones(T,n_batch,n_ctrl),
                        lqr_iter=20,
                        verbose=1,
                        backprop=False,
                        exit_unconverged=False,
                        grad_method=mpc.GradMethods.AUTO_DIFF)
    
    xs, us, objs = ctrl_mpc(x_init, QuadCost(C,c), network)
    print(xs)
    