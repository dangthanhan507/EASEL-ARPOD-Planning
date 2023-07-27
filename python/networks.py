import torch.nn as nn
import torch
from arpod_dynamics import createHCWDiscreteMatrices

# class ModifiedControlNetwork(nn.Module):
#     def __init__(self, T, mu_GM, R):
#         super().__init__()
#         self.A, self.B = createHCWDiscreteMatrices(T,mu_GM,R)

#         self.layers = []
#         self.layers.append(nn.Linear(9,64))
#         self.layers.append(nn.Linear(64,128))
#         self.layers.append(nn.Linear(128,256))
#         self.layers.append(nn.Linear(256,256))
#         self.layers.append(nn.Linear(256,256))
#         self.layers.append(nn.Linear(256,128))
#         self.layers.append(nn.Linear(128,64))
#         self.layers.append(nn.Linear(64,8))
#         self.model = nn.Sequential(*self.layers)

#     def forward(self, state, action):
#         pass

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
    
    def forward(self, inp):
        '''
        '''
        return self.model(inp)