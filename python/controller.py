#=========================================
# controller.py
# 
# 
#=========================================
from arpod_dynamics import createHCWMatrices, createHCWDiscreteMatrices
import control
import numpy as np

import casadi   as cs
import l4casadi as l4c

'''
    LQR utilized for ARPOD rendezvous
'''
class LQR:
    def __init__(self, A, B, Q, R):
        '''
        '''
        self.K,self.S,self.E = control.lqr(A,B,Q,R)
    def getControl(self, x):
        '''
            x = (6,1)
        '''

        # (3,6) @ (6,1) = (3,1) which is control input
        u = (-self.K) @ x

        #set limits
        u[u>=0.1] = 0.1
        u[u<=-0.1] = -0.1
        return u

'''
NOTES:
    SX is symbolic expression for Symbolic differentiation (limited to binary, unary, and scalar operation)
    MX is matrix expression (not limited)
'''
class NN_MPC:
    #MPC using neural network as complete replacement for dynamics
    #max thrust is 0.1
    def __init__(self, tstep, Q, R, horizon, network):
        self.dt = tstep
        self.horizon = horizon
        self.network = network
        self.Q = Q
        self.R = R

        self.x_dim = Q.shape[0]
        self.u_dim = R.shape[0]

    def fn_objective(self, x, u, Q, R):
        #discrete objective
        L = cs.sum1(cs.SX(np.diag(Q))*(x*x)) + cs.sum1(cs.SX(np.diag(R))*(u*u))
        return L
    def fn_constr(self):
        pass

    def extract_solution(self,sol):
        solution  = np.array(sol)
        final_x_idx = self.x_dim*self.horizon

        xs = solution[0:final_x_idx]
        us = solution[final_x_idx:]
        x = xs.reshape((6,-1))
        u = us.reshape((3,-1))
        return x,u
    def optimize(self, state):
        u = cs.SX.sym('u',self.u_dim,self.horizon)
        x = cs.SX.sym('x',self.x_dim,self.horizon)
        # x0 = cs.SX(state)
        L = self.fn_objective(x,u,self.Q,self.R)

        #setup optimization problem (objective + constraints) + solver

        opt_x  = cs.vertcat(x.reshape((-1,1)),u.reshape((-1,1)))
        nlp = cs.nlpsol('nlp','ipopt', {'x': opt_x, 'g': })

        #solve
        solution = nlp()
        x,u = self.extract_solution(solution['x'])
        return x,u




if __name__ == '__main__':
    x0 = np.array([[-10,10,10,1,1,1]]).T
    mu_GM = 100 
    R = 1000
    tstep = 1
    A,B = createHCWMatrices(tstep,mu_GM,R)
    Ak,Bk = createHCWDiscreteMatrices(tstep, mu_GM, R)

    Q = np.eye(6)
    R = np.eye(3)

    controller = LQR(A,B,Q,R)

    x = x0.copy()
    for iter in range(1000):
        u = controller.getControl(x)
        x = Ak@x + Bk@(u)

    print(x)

