import numpy as np

def get_mse(A,X,ntimepts,nreps):
    X = X[:,:ntimepts].reshape(len(X),(ntimepts)*nreps,order='F')
    X_pred = np.zeros((A.shape[0],ntimepts*nreps))
    count = 0
    for i in range(0,nreps):
        x_test_ic = X[:,i*(ntimepts):i*(ntimepts)+1]
        for j in range(0,ntimepts):
            X_pred[:,count:count+1] = np.dot(np.linalg.matrix_power(A,j),x_test_ic) 
            count += 1
    mse = 1/(ntimepts*nreps) * np.linalg.norm(X_pred - X,ord=2)**2 / np.linalg.norm(X,ord=2)
    return mse

def n_step_prediction(A,X,ntimepts,nreps):
    X = X[:,:ntimepts].reshape(len(X),(ntimepts)*nreps,order='F')
    X_pred = np.zeros((A.shape[0],ntimepts*nreps))
    count = 0
    for i in range(0,nreps):
        x_test_ic = X[:,i*(ntimepts):i*(ntimepts)+1]
        for j in range(0,ntimepts):
            X_pred[:,count:count+1] = np.dot(np.linalg.matrix_power(A,j),x_test_ic) 
            count += 1
    feature_means = np.mean(X,axis=1).reshape(len(X),1)
    cd = 1 - ((np.linalg.norm(X - X_pred,ord=2)**2)/(np.linalg.norm(X - feature_means,ord=2)**2))   # coeff of determination aka R^2 
    print(f'Coefficient of determination for n-step prediction is {cd:.3e}')
    return X_pred, cd

def trim_weights(A,thresh):
    nnonzero = len(np.nonzero(A)[0])
    Atrim = (np.absolute(A) > thresh) * A
    nnonzero_sparse = len(np.nonzero(Atrim)[0])
    percent_nonzero_to_zero = (nnonzero - nnonzero_sparse)/nnonzero
    # print('Percent of nonzero parameters set to zero:',percent_nonzero_to_zero)
    return Atrim, percent_nonzero_to_zero


def dmd(data,rank_reduce=True,r=None,trim=False,trimThresh=1.5e-3):
    '''
    data (shape nxmxr): features along rows, time along columns, trajectories along depth 
    if trim: 
        trimThresh: any element of A that is less than trimThresh will be set to zero
    '''

    # form the snapshot matrices for dmd
    Xp = data[:,:-1].reshape(len(data),(data.shape[1]-1)*data.shape[2],order='F') 
    Xf = data[:,1:].reshape(len(data),(data.shape[1]-1)*data.shape[2],order='F')
#     # to add bias term, uncomment below
#     ones_row = np.ones((1,Xp.shape[1]))
#     Xp = np.concatenate((Xp,ones_row),axis=0)
#     Xf = np.concatenate((Xf,ones_row),axis=0)
#     ones_row = np.ones((1,Xp.shape[1]+2))
#     nR = data.shape[2] # number of replicates
#     data = np.concatenate((data.reshape(len(data),data.shape[1]*data.shape[2],order='F'),ones_row),axis=0)
#     data = data.reshape(len(data),data.shape[1]//nR,nR,order='F')

    # exact dmd algorithm
    if rank_reduce == False:
        A = Xf @ np.linalg.pinv(Xp)
    else: 
        U,s,Vh = np.linalg.svd(Xp)
        if rank_reduce and r==None: 
            r = np.minimum(U.shape[1],Vh.shape[0])
        U_r = U[:,0:r] # truncate to rank-r
        s_r = s[0:r]
        Vh_r = Vh[0:r,:]
        Atilde = U_r.T @ Xf @ Vh_r.T @ np.diag(1/s_r) # low-rank dynamics
        A = U_r@Atilde@U_r.T # full model A

    # make model sparse 
    if trim:
        A, percent_nonzero_to_zero = trim_weights(A,trimThresh)

    # calculate prediction accuracy 
#     X_pred, cd = n_step_prediction(A,data,data.shape[1],data.shape[2])
#     X_pred = X_pred.reshape(len(data),data.shape[1],data.shape[2])
    data_red = np.zeros((r,data.shape[1],data.shape[2]))
    data_red[:,:,0] = np.dot(U_r.T ,data[:,:,0])
    data_red[:,:,1] = np.dot(U_r.T ,data[:,:,1])
    X_pred_red, cd_red = n_step_prediction(Atilde,data_red,data_red.shape[1],data_red.shape[2])
    X_pred_red = X_pred_red.reshape(len(data_red),data_red.shape[1],data_red.shape[2])
    
    # compute eigenvectors and eigenvalues of DMD operator
    L,W = np.linalg.eig(Atilde)
    
    # compute DMD modes
    Phi = Xf @ Vh_r.T @ np.diag(1/s_r) @ W
    
    # compute mode amplitudes for the two replicates (or temporal modes)
#     b_r0 = np.linalg.pinv(Phi) @ data[:,0:1,0]
#     b_r1 = np.linalg.pinv(Phi) @ data[:,0:1,1]
#     b_r0 = np.linalg.pinv(Phi) @ X_pred[:,:,0]
#     b_r1 = np.linalg.pinv(Phi) @ X_pred[:,:,1]
    b_r0 = np.linalg.inv(np.dot(W,np.diag(L))) @ data_red[:,:,0]
    b_r1 = np.linalg.inv(np.dot(W,np.diag(L))) @ data_red[:,:,1]

    return A, Atilde, cd_red, L, W, Phi, b_r0, b_r1









