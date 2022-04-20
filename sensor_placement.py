import numpy as np
from scipy.optimize import minimize
from scipy import stats 
import time
from copy import deepcopy

def objective_manyIC(x,At,x0):
    '''
    The objective is to maximize the output energy
    x is the design variable or the output gain
    At is the KO raised to the power t
    x0 is the initial condition
    '''
    obj = 0
    for i in range(len(At)):
        obj += -np.linalg.norm(x@At[i]@x0,ord=2) # for use with energy_max or energy_max_uni_dist
#         obj += -np.linalg.norm( (x@At[i]@x0).reshape(1,1) ,ord=2) # for use with energy_max_with_mean
    return obj

def energy_maximization_with_uni_dist(X,A,ntimepts,Tf,IC=0):    
    print('------Optimizing for gene sampling weights------')
    start_time = time.time()
    r = np.random.randint(0,100)
    np.random.seed(r)
    C0 = np.random.normal(4000,1000.0,size=(1,len(A))) # initialization
    At = []
    for i in range(0,Tf+1):
        At.append(np.linalg.matrix_power(A,i))
    x0 = X[:,0,:]
    # get the min and max of each gene's IC
    x0min = np.min(x0,axis=1) # x0min = np.min(X[:,0,:],axis=1)
    x0max = np.max(x0,axis=1) # x0max = np.max(X[:,0,:],axis=1)
    # form a set of IC's distributed uniformly
    numICs = X.shape[0]
    x0uni = np.zeros((len(x0min),numICs))
    for ii in range(x0uni.shape[1]):
        x0tmp = np.random.uniform(x0min,x0max)
        x0uni[:,ii] = x0tmp
    print('Initial objective: ' + str(objective_manyIC(C0,At,x0uni)))
    # optimize
    bnds = tuple([(0.0,None) for i in range(X.shape[0])]) # C should be nonnegative
    solution = minimize(objective_manyIC,C0,args=(At,x0uni),method='SLSQP',bounds=bnds)
    C = (solution.x).reshape(1,C0.shape[1])
    # show final objective
    print('Final objective: ' + str(objective_manyIC(C,At,x0uni)))
#     C = (C/C.max()).T # normalizing C to be b/w 0 and 1
    print((time.time() - start_time)/60, 'minutes')

    return C,r

def reconstruct_x0(data_fc_norm,nT,A,C,CsortedInds,samplingFreq=5,order='top'):

    '''
    nT is the number of timepoints for which to generate outputs and reconstruct x0 with.
    First sample genes by rank (or randomly), compute reconstruct error, 
    then sample again (adding to the already sampled genes)
    and repeat
    '''
    
    topNums = list(range(0,len(data_fc_norm),samplingFreq)) # sampling these genes in order of rank samplingFreq at a time
    # don't need to sample all genes, as we know the reconstruction loss asymptotes after 80 genes or so.
    rho = []
    Uset = list(range(len(data_fc_norm))) # set of all genes
    NUset = deepcopy(Uset) # updated each iteration to remove the genes that were used prior (no double dipping)
    for ii in range(len(topNums)-1): 
        if order == 'top':
            top_inds = list(np.array(CsortedInds[::-1])[topNums[0]:topNums[ii+1]])
            top_inds_complement = list(set(list(range(len(data_fc_norm)))) - set(top_inds))
        elif order == 'random':
            if len(NUset) <= samplingFreq: 
                break 
            if ii == 0: 
                rs = random.sample(NUset,samplingFreq)
                store_rs = deepcopy(rs)
            else: 
                rs = random.sample(NUset,samplingFreq)
                store_rs = store_rs + rs
            NUset = list(set(NUset) - set(store_rs))
            top_inds_complement = deepcopy(NUset)
            
        C_sensors = deepcopy(np.array(C).reshape(1,-1))
        C_sensors[:,top_inds_complement] = 0.0

        # generate outputs at nT timepoints, check rank of O_T and estimate IC
        O_T = np.zeros((len(C_sensors)*nT,C_sensors.shape[1])) # O_T has dim nOutputs*nTimepoints,nStateVars
        for ii in range(O_T.shape[0]): # observability matrix
            O_T[ii] = C_sensors @ np.linalg.matrix_power(A,ii)
#         print('Rank of observability matrix: ',np.linalg.matrix_rank(O_T)) # don't print, too slow

        # use the learned dynamics and sampling to generate output over time
        y_r1 = np.zeros((1,nT)) 
        y_r2 = np.zeros((1,nT))
        for ii in range(nT):
            y_r1[:,ii] = C_sensors @ np.linalg.matrix_power(A,ii) @ data_fc_norm[:,0,0]
            y_r2[:,ii] = C_sensors @ np.linalg.matrix_power(A,ii) @ data_fc_norm[:,0,1]

        # estimate of x0 
        x0_est_r1 = np.linalg.pinv(O_T) @ y_r1.T
        x0_est_r2 = np.linalg.pinv(O_T) @ y_r2.T

        rho1 = stats.pearsonr(data_fc_norm[:,0,0],np.squeeze(x0_est_r1))[0]
        rho2 = stats.pearsonr(data_fc_norm[:,0,1],np.squeeze(x0_est_r2))[0]
        rho.append((rho1 + rho2)/2)
        
    return rho

def gram_matrix(A,x0,nT=50,reduced=True,projection_matrix=np.array([])):

    '''
    A: matrix representation of the Koopman operator
    x0: initial conditions from measurements
    nT: number of timepoints over which to compute the Gram matrix
    reduced: if True, will compute reduced G from reduced data and KO and will also return full G after inverse projection
    projection_matrix: the matrix used to data and KO to low-dimensional space (first r eigenvectors of Data.T @ Data)
    Both A and x0 can be either the full dimensional data and KO or they can be the DMD projected data and KO
    If projected, then return both the projected G and the full G after inverting the projection
    If not projected, then compute full G (can be slow, especially if the data dimension exceeds a couple thousand)
    Furthermore, for sensor placement we need to compute the eigendecomposition of G, so having the reduced G is handy   
    '''

    # generate artificial initial conditions for robust optimization 
    # get the min and max of each gene's initial value
    x0min = np.min(x0,axis=1)
    x0max = np.max(x0,axis=1)
    # form a set of new initial conditions distributed uniformly from x0min to x0max
    numICs = x0.shape[0]
    x0uni = np.zeros((len(x0min),numICs))
    x0uni[:,0:x0.shape[1]] = deepcopy(x0)
    for ii in range(x0.shape[1],x0uni.shape[1]):
        x0tmp = np.random.uniform(x0min,x0max)
        x0uni[:,ii] = x0tmp
    G = np.zeros_like(A)
    for ii in range(nT):
        A_pow = np.linalg.matrix_power(A,ii)
        G += np.matmul( np.matmul(A_pow, x0uni), np.matmul(x0uni.T, A_pow.T) ) 
    # right eigenvectors of G (columns of V) are rows of the gene sampling matrix (or vector if just one eigvec kept)

    if reduced: 
        Gfull = np.matmul(np.matmul(projection_matrix, G), projection_matrix.T)
        return G, Gfull
    else: 
        return G # this is the full G, computed directly from full KO and data












