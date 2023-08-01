#=========================================
# controller.py
# 
# 
#=========================================
from arpod_dynamics import createHCWMatrices, createHCWDiscreteMatrices
import control
import numpy as np

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

