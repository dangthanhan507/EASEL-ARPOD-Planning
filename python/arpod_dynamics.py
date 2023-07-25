import torch.nn as nn
import torch
import numpy as np

'''
    Research references:
    --------------------

        MAML:
        https://arxiv.org/pdf/1703.03400.pdf
        https://github.com/cbfinn/maml

        Adapt MBRL
        https://arxiv.org/pdf/1803.11347.pdf
        https://github.com/iclavera/learning_to_adapt/tree/master

'''
def createHCWMatrices(T, mu_GM, R):
    n = np.sqrt(mu_GM / (R*R*R))
    A = np.zeros(6,6)
    B = np.zeros(6,3)

    S = np.sin(n*T)
    C = np.cos(n*T)
    
    A[0,:] = np.array([4-3*C,0,0,S/n,2*(1-C)/n,0])
    A[1,:] = np.array([6*(S-n*T),1,0,-2*(1-C)/n,(4*S-3*n*T)/n,0])
    A[2,:] = np.array([0,0,C,0,0,S/n])
    A[3,:] = np.array([3*n*S,0,0,C,2*S,0])
    A[4,:] = np.array([-6*n*(1-C),0,0,-2*S,4*C-3,0])
    A[5,:] = np.array([0,0,-n*S,0,0,C])

    B[0,:] = np.array([(1-C)/(n*n),(2*n*T-2*S)/(n*n),0])
    B[1,:] = np.array([-(2*n*T-2*S)/(n*n),(4*(1-C)/(n*n))-(3*T*T/2),0])
    B[2,:] = np.array([0,0,(1-C)/(n*n)])
    B[3,:] = np.array([S/n,2*(1-C)/n, 0])
    B[4,:] = np.array([-2*(1-C)/n,(4*S/n) - (3*T),0])
    B[5,:] = np.array([0,0,S/n])

    return A,B


class ModifiedControlNetwork(nn.Module):
    def __init__(self, T, mu_GM, R):
        super().__init__()
        self.A, self.B = createHCWMatrices(T,mu_GM,R)

        self.layers = []
        self.layers.append(nn.Linear(9,64))
        self.layers.append(nn.Linear(64,128))
        self.layers.append(nn.Linear(128,256))
        self.layers.append(nn.Linear(256,256))
        self.layers.append(nn.Linear(256,256))
        self.layers.append(nn.Linear(256,128))
        self.layers.append(nn.Linear(128,64))
        self.layers.append(nn.Linear(64,8))
        self.model = nn.Sequential(*self.layers)

    def forward(self, state, action):
        pass

class DynamicsNetwork(nn.Module):
    def __init__(self):
        super().__init__() #must call

        #setting up hidden layer

        # control input is 3, state input is 6... that makes 9
        self.layers = []
        self.layers.append(nn.Linear(9,64))
        self.layers.append(nn.Linear(64,128))
        self.layers.append(nn.Linear(128,256))
        self.layers.append(nn.Linear(256,256))
        self.layers.append(nn.Linear(256,256))
        self.layers.append(nn.Linear(256,128))
        self.layers.append(nn.Linear(128,64))
        self.layers.append(nn.Linear(64,8))
        self.model = nn.Sequential(*self.layers)
    
    def forward(self, state, action):
        '''
        '''
        inp = torch.cat(state,action)
        return self.model(inp)
    

