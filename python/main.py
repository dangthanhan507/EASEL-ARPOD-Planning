from arpod_dynamics import createHCWDiscreteMatrices
from controller import MPC
import numpy as np
if __name__ == '__main__':
    A,B = createHCWDiscreteMatrices(1,1000,100)
    D = np.array([[1, -1,  0,  0,  0,  0],
                  [0,  0,  1, -1,  0,  0],
                  [0,  0,  0,  0,  1, -1]])
    Q = np.eye(6)*100
    R = np.eye(6)

    pos0 = np.ones(3)*-5
    vel0 = np.ones(3)*0.1
    x0 = np.array([pos0,vel0]).reshape((6,1))
    mpc = MPC(1,Q,R,horizon=5)
    
    t = 0
    finalT = 100
    while t < finalT:
        xs,us = mpc.optimize(x0)
        u = us[:,0].reshape((6,1))
        x0 = A@x0 + B@D@u

        t += 1

    print(x0)