import sys
sys.path.append('./ml-casadi/')
import casadi   as cs
import torch
import torch.nn as nn
import numpy as np

import l4casadi as l4c
import ml_casadi.torch as mc

def createGTMatrices():
    '''
    '''
    T = 1
    mu_GM = 1000
    R = 100

    n = np.sqrt(mu_GM / (R*R*R))
    A = np.zeros((6,6))
    B = np.zeros((6,3))

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
        return self.model(inp.T).T

if __name__ == '__main__':
    device = 'cpu'
    network = SimpleDynamicsNetwork()
    network_load = torch.load('./models/hcw_simple_model.pt').to(device)
    network.load_state_dict(network_load.state_dict())

    ml_model = mc.TorchMLCasadiModuleWrapper(network,input_size=9,output_size=6)

    xu = cs.MX.sym('xu', 9)
    output = ml_model.approx(xu, order=1)
    cs_func = cs.Function('model_approx_wrapper', [xu, ml_model.sym_approx_params(order=1,flat=True)], [output])

    Q = np.eye(6)
    R = np.eye(3)

    pos0 = np.ones(3)*-5
    vel0 = np.ones(3)*0.1
    x0_np = np.array([pos0,vel0]).reshape((6,1))
    x0 = cs.MX(x0_np)
    x = cs.MX.sym('x',6,10)
    u = cs.MX.sym('u',3,10)

    #setup objective
    L = cs.sum2(cs.sum1(cs.MX(np.diag(Q))*(x*x)) + cs.sum1(cs.MX(np.diag(R))*(u*u)))

    #setup dynamic constraint
    xs = cs.horzcat(x0,x[:,:-1])

    #setup using neural net
    outs_net = []
    for i in range(xs.shape[1]):
        inp_net = cs.vertcat(xs[:,i], u[:,i])

        casadi_param = ml_model.approx_params(inp_net, order=1, flat=True)
        out = cs_func(inp_net,casadi_param)
        outs_net.append(out)
    pred_xs = cs.horzcat(*outs_net)

    #setup using ground truth dynamics
    # pred_xs = A@xs + B@u

    dyn_constr = (pred_xs - x).reshape((-1,1))

    #setup decision variable constraints
    m,n = x.shape
    ubx_x = [cs.inf]*(m*n)
    lbx_x = [-cs.inf]*(m*n)
    m,n = u.shape
    ubx_u = [0.1]*(m*n)
    lbx_u = [-0.1]*(m*n)

    ubx = ubx_x + ubx_u
    lbx = lbx_x + lbx_u

    #setup decision variables
    opt_x = cs.vertcat(x.reshape((-1,1)), u.reshape((-1,1)))

    #solve
    options = {'printLevel': "none"}
    solver = cs.qpsol('solver','qpoases', {'f': L, 'x': opt_x, 'g': dyn_constr}, options)
    solution = solver(lbx=lbx,ubx=ubx,lbg=0,ubg=0)
    print('\n\ndynamics constraint output:', solution['g']) 


    #use solution to get optimized control input and state vars
    decision_solution  = np.array(solution['x'])
    final_x_idx = 6*10 # dim * horizon
    x_sol = decision_solution[:final_x_idx].reshape((-1,6)).T
    u_sol = decision_solution[final_x_idx:].reshape((-1,3)).T

    print('\n\sol output1:')
    print(x_sol[:,0])

    print('l4casadi output:')
    out = l4c_model(np.vstack((x0_np,u_sol[:,0].reshape((3,1)))))
    print(out)

    net_x = torch.Tensor(np.vstack((x0_np,u_sol[:,0].reshape((3,1)))))
    out = network(net_x)
    print("torch output1:")
    print(out)

    print('sol output2:')
    print(x_sol[:,1])

    net_x = torch.Tensor(np.vstack((x_sol[:,0].reshape((6,1)),u_sol[:,1].reshape((3,1)))))
    out = network(net_x)
    print('torch output2:')
    print(out)

    print('sol output3:')
    print(x_sol[:,2])

    net_x = torch.Tensor(np.vstack((x_sol[:,1].reshape((6,1)),u_sol[:,2].reshape((3,1)))))
    out = network(net_x)
    print('torch output3:')
    print(out)