import numpy as np
import datatable as dt
import pandas as pd
from scipy.sparse import issparse, csr_matrix, find, coo_matrix
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.decomposition import PCA
from scipy.spatial.distance import squareform
from sklearn.neighbors import NearestNeighbors
from scipy.sparse.linalg import eigs
from numpy.linalg import norm


def magicBatch(data=None, mar_mat_input=None, n_pca_components=20, n_diffusion_components = 0, random_pca=True,t=None, 
          k=30, ka=10, epsilon=1, rescale_percent = 90, rescale_method = 'classic', csv_l=None, csv_d=None, csv_i=None, to_print=None):
    if(to_print != None):
        print(to_print)
    if(isinstance(data, str)):
        data = pd.read_csv(data)
        if data.iloc[0,0] == 'fzJKy2SHmdwIPQ':
            markov_mat_input = run_pca(data.iloc[1:, np.where(data.iloc[0,:] == 1)[0]], n_components=n_pca_components, random=random_pca)
            data = data.iloc[1:,1:]
    if not isinstance(mar_mat_input, str) and not isinstance(mar_mat_input, np.ndarray):
        print('doing PCA')
        markov_mat_input = run_pca(data, n_components=min(n_pca_components, data.shape[0]), random=random_pca)
    else:
        if not isinstance(mar_mat_input, np.ndarray):
            print('using user-supplied dimensional reduction')
            markov_mat_input = mar_mat_in = pd.read_csv(mar_mat_input)
        else:
            markov_mat_input = mar_mat_input

    L = compute_markov(markov_mat_input, k=k, epsilon=epsilon, 
                       distance_metric='euclidean', ka=ka)
    if csv_l != None:
        L_coo = coo_matrix(L)
        L_coo = dt.Frame(row=L_coo.row, col=L_coo.col, data=L_coo.data)
        L_coo.to_csv(csv_l)

    if n_diffusion_components > 0:                
        # Eigen value decomposition
        D, V = eigs(L, k=n_diffusion_components, tol=1e-4, maxiter=1000)
        D = np.real(D)
        V = np.real(V)
        inds = np.argsort(D)[::-1]
        D = D[inds]
        V = V[:, inds]
        V = L.dot(V)
        for i in range(V.shape[1]):
            V[:, i] = V[:, i] / norm(V[:, i])
        V = np.round(V, 10)
        V = dt.Frame(V)
        V.to_csv(csv_d)

    if t != None:
        L = L.todense()
        if(isinstance(t, str)):
            t_param = t.split('_')
            t_param = [int(i) for i in t_param]
        data_new = multi_t_fast_impute(data, L = L, t = t_param, rescale_percent = rescale_percent, rescale_method = rescale_method)
        
    if csv_i != None:
        data_new = dt.Frame(np.hstack(data_new))
        data_new.to_csv(csv_i)

def run_pca(data, n_components=100, random=True):

    solver = 'randomized'
    if random != True:
        solver = 'full'

    pca = PCA(n_components=n_components, svd_solver=solver)
    return pca.fit_transform(data)

def multi_t_fast_impute(data, L, t, rescale_percent = 0, rescale_method = 'adaptive'):

    t_max = max(t)
	
    if t_max >= 2: 
        L2 = np.matmul(L, L)
    if t_max >= 4: 
        L4 = np.matmul(L2, L2)
    if t_max >= 8: 
        L8 = np.matmul(L4, L4)
    if t_max >= 16: 
        L16 = np.matmul(L8, L8)
    if t_max >= 32:
        L32 = np.matmul(L16, L16)
    if len(set([3,7,11,15,19,23,27,31,35]) & set(t)) > 0:
        L3 = np.matmul(L2, L)

    dif_op = ['empty',
              'L',
              'L2',
              'L3',
              'L4',
              'np.matmul(L4, L)',
              'np.matmul(L4, L2)',
              'np.matmul(L4, L3)',
              'L8',
              'np.matmul(L8, L)',
              'np.matmul(L8, L2)',
              'np.matmul(L8, L3)',
              'np.matmul(L8, L4)',
              'np.matmul(np.matmul(L8, L4), L)',
              'np.matmul(np.matmul(L8, L4), L2)',
              'np.matmul(np.matmul(L8, L4), L3)',
              'L16',
              'np.matmul(L16, L)',
              'np.matmul(L16, L2)',
              'np.matmul(L16, L3)',
              'np.matmul(L16, L4)',
              'np.matmul(np.matmul(L16, L4), L)',
              'np.matmul(np.matmul(L16, L4), L2)',
              'np.matmul(np.matmul(L16, L4), L3)',
              'np.matmul(L16, L8)',
              'np.matmul(np.matmul(L16, L8), L)',
              'np.matmul(np.matmul(L16, L8), L2)',
              'np.matmul(np.matmul(L16, L8), L3)',
              'np.matmul(np.matmul(L16, L8), L4)',
              'np.matmul(np.matmul(np.matmul(L16, L8), L4), L)',
              'np.matmul(np.matmul(np.matmul(L16, L8), L4), L2)',
              'np.matmul(np.matmul(np.matmul(L16, L8), L4), L3)',
              'L32',
              'np.matmul(L32, L)',
              'np.matmul(L32, L2)',
              'np.matmul(L32, L3)',
              'np.matmul(L32, L4)']
    l_vars = locals()
    L_t = [eval(dif_op[i], {'np': np}, l_vars) for i in t]
    print('MAGIC: data_new = L_t * data')
    data_new = [np.array(np.dot(i, data)) for i in L_t]
    if rescale_percent > 0:
        data_new = [rescale_data(data, imputed_data = i, rescale_percent = rescale_percent, rescale_method = rescale_method) for i in data_new]
    return data_new
    
    
def rescale_data(data, imputed_data, rescale_percent, rescale_method):
    print('Rescaling Data')
    rescale_percent = rescale_percent / 100
    num_rows = data.shape[0]
    num_cols = data.shape[1]
    if rescale_method == "classic":
        print('using classic method')
        M99 = np.percentile(data, rescale_percent, axis=0)
        M100 = data.max(axis=0)
        indices = np.where(M99 == 0)[0]
        M99[indices] = M100[indices]
        M99_new = np.percentile(imputed_data, rescale_percent, axis=0)
        M100_new = imputed_data.max(axis=0)
        indices = np.where(M99_new == 0)[0]
        M99_new[indices] = M100_new[indices]
    if rescale_method == "adaptive":
        print('using adaptive method')
        nz = np.count_nonzero(data, 0) / num_rows
        z = 1 - nz
        rescale_percent_adj = rescale_percent * nz + z
        M99 = np.asarray([np.quantile(data.iloc[:,i], rescale_percent_adj[i], axis=0) for i in range(num_cols)])
        M100 = data.max(axis=0)
        indices = np.where(M99 == 0)[0]
        M99[indices] = M100[indices]
        M99_new = np.asarray([np.quantile(imputed_data[:,i], rescale_percent_adj[i], axis=0) for i in range(num_cols)])
        M100_new = imputed_data.max(axis=0)
        indices = np.where(M99_new == 0)[0]
        M99_new[indices] = M100_new[indices]
    np.seterr(invalid='ignore')
    max_ratio = np.divide(M99, M99_new)
    indices = np.isnan(max_ratio)
    max_ratio[indices] = 0
    imputed_data = np.multiply(imputed_data, np.tile(max_ratio, (len(data), 1)))
    return imputed_data


def compute_markov(data, k=10, epsilon=1, distance_metric='euclidean', ka=0):
    # This function is identical the "compute_markov" function from the classic implementation:
    # https://github.com/dpeerlab/magic
    N = data.shape[0]
    k = min(k, data.shape[0] -1)
    # Nearest neighbors
    print('Computing distances')
    nbrs = NearestNeighbors(n_neighbors=k, metric=distance_metric).fit(data)
    distances, indices = nbrs.kneighbors(data)

    if ka > 0:
        print('Autotuning distances')
        for j in reversed(range(N)):
            temp = sorted(distances[j])
            lMaxTempIdxs = min(ka, len(temp))
            if lMaxTempIdxs == 0 or temp[lMaxTempIdxs] == 0:
                distances[j] = 0
            else:
                distances[j] = np.divide(distances[j], temp[lMaxTempIdxs])

    # Adjacency matrix
    print('Computing kernel')
    rows = np.zeros(N * k, dtype=np.int32)
    cols = np.zeros(N * k, dtype=np.int32)
    dists = np.zeros(N * k)
    location = 0
    for i in range(N):
        inds = range(location, location + k)
        rows[inds] = indices[i, :]
        cols[inds] = i
        dists[inds] = distances[i, :]
        location += k
    if epsilon > 0:
        W = csr_matrix( (dists, (rows, cols)), shape=[N, N] )
    else:
        W = csr_matrix( (np.ones(dists.shape), (rows, cols)), shape=[N, N] )

    # Symmetrize W
    W = W + W.T

    if epsilon > 0:
        # Convert to affinity (with selfloops)
        rows, cols, dists = find(W)
        rows = np.append(rows, range(N))
        cols = np.append(cols, range(N))
        dists = np.append(dists/(epsilon ** 2), np.zeros(N))
        W = csr_matrix( (np.exp(-dists), (rows, cols)), shape=[N, N] )

    # Create D
    D = np.ravel(W.sum(axis = 1))
    D[D!=0] = 1/D[D!=0]

    #markov normalization
    T = csr_matrix((D, (range(N), range(N))), shape=[N, N]).dot(W)

    return T


