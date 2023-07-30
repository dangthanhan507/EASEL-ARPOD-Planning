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

    neural networks only use MX so we stuck with that :(
    
    TODO:
        make it so that I can swap out the dynamic constraints easily
'''
class NN_MPC:
    #MPC using neural network as complete replacement for dynamics
    #max thrust is 0.
    def __init__(self, tstep, Q, R, horizon, network):
        self.dt = tstep
        self.horizon = horizon
        self.network = l4c.L4CasADi(network,has_batch=False,device='cpu')
        self.Q = Q
        self.R = R

        self.x_dim = Q.shape[0]
        self.u_dim = R.shape[0]

    def fn_objective(self, Q, R):
        #convert x_fn into usable function for objective

        x = cs.MX.sym('X', self.x_dim, self.horizon)
        u = cs.MX.sym('U', self.u_dim, self.horizon)

        J = cs.sum2(cs.sum1(cs.MX(np.diag(Q))*(x*x)) + cs.sum1(cs.MX(np.diag(R))*(u*u)))
        L = cs.Function('loss_fn', [x,u], [J])
        return L #return loss function
    def prop_network(self, x, u):

        #stack x,u
        inp_net = cs.vertcat(x, u)
        out = self.network(inp_net.T)

        #(6,N)
        return out.T
    
    def fn_dynamic_constr(self, x0, x, u):
        xs = cs.horzcat(x0,x[:,:-1]) # (x0.... xn-1)

        dyn_x = self.prop_network(xs,u) # predicted (x1... xn)

        # xk+1 = f(xk,uk)
        # returns f(xk,uk)-xk+1 = 0 flattened 
        return (x - dyn_x).reshape((-1,1))
    
    def fn_decision_constr(self, x, u):
        #setup decision variable limits
        m,n = x.shape
        ubx = (np.ones(m*n)*np.inf).tolist()
        lbx = (np.ones(m*n)*-np.inf).tolist()
        m,n = u.shape
        ubu = (np.ones(m*n)*0.1).tolist()
        lbu = (np.ones(m*n)*-0.1).tolist()

        ub_d = ubx + ubu
        lb_d = lbx + lbu
        return ub_d, lb_d

    def extract_solution(self,sol):
        solution  = np.array(sol)
        final_x_idx = self.x_dim*self.horizon

        xs = solution[0:final_x_idx]
        us = solution[final_x_idx:]
        x = xs.reshape((6,-1))
        u = us.reshape((3,-1))
        return x,u
    def optimize(self, state):
        u = cs.MX.sym('u',self.u_dim,self.horizon)
        x = cs.MX.sym('x',self.x_dim,self.horizon)
        x0 = cs.MX(state)

        #setup optimization problem (objective + constraints) + solver
        '''
            NOTE: General Structure of nlpsol

            min f(x)
            s.t. lbg <= g(x) <= ubg
                 lbx <=   x  <= ubx


            for our purposes, we should do some work on our end to prevent
            extra complications modifying lbg and ubg. Normally, these params
            are meant for how far we can adjust our decision variables.

            for equality constraints like:
                xk+1 = A*xk + B*uk
            lets leave it as
                0 = A*xk + B*uk - xk+1

            so we can let: lb, ub = 0
        '''
        L = self.fn_objective(self.Q,self.R)
        dyn_x  = self.fn_dynamic_constr(x0,x,u)
        # lbg    = cs.MX(np.zeros(dyn_x.shape[0]))
        # ubg    = cs.MX(np.zeros(dyn_x.shape[0]))
        ubx,lbx = self.fn_decision_constr(x,u)
        opt_x  = cs.vertcat(x.reshape((-1,1)),u.reshape((-1,1)))    

        options = {'ipopt.print_level': 1}
        solver = cs.nlpsol('solver','ipopt', {'f': L(x,u), 'x': opt_x, 'g': dyn_x}, options)
        print(solver)

        initial_guess = cs.MX(np.zeros(opt_x.shape[0]))
        #solve
        solution = solver(lbx=lbx,ubx=ubx,lbg=0,ubg=0)
        
        print(solution['x'])
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

