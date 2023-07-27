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
def createHCWDiscreteMatrices(T, mu_GM, R):
    '''
    '''
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

def createHCWMatrices(T, mu_GM, R):
    '''
        Create HCW matrices
    '''
    n = np.sqrt(mu_GM / (R*R*R))


    #creating A matrix
    A = np.hstack((np.zeros((3,3)),np.eye(3)))
    #bottom left and bottom right part of A matrix
    BL = np.array([[3*n*n, 0, 0],[0,0,0],[0,0,-n*n]])
    BR = np.array([[0,2*n,0],[-2*n,0,0],[0,0,0]])
    lowerA = np.hstack((BL,BR))
    A = np.vstack((A,lowerA))
    B = np.vstack((np.zeros((3,3)), np.eye(3)))

    return A,B