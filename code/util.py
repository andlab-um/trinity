# -*- coding: utf-8 -*-

__author__ = 'Zhaoning Li'

import warnings
import pandas as pd
import pickle
import numpy as np
from scipy.stats.mstats_basic import rankdata
from scipy.spatial import distance
from tqdm import tqdm
import scipy.sparse as sp
from scipy.optimize import linprog

warnings.filterwarnings('ignore')
exclude_num = [103, 112, 113, 115, 125, 128, 131]


def excluded_pair(order):

    """
    Find non-informative subject pairs (comparisons of subjects to themselves) to be exculded

    Parameters
    ----------
    order: array.
        Number of the subject.

    Returns
    -------
    result: list.
        The list of tuples, i.e., non-informative subject pairs to be exculded.
    """
    
    result = []
    for i, item_x in enumerate(order):
        for j, item_y in enumerate(order):
            if i < j and item_x == item_y:
                result.append((i, j))
                
    return result


def wasserstein_distance(p, q, D):

    """
    Calculate Wasserstein distance by linear programming
    Adapted from: Jianlin Su @ https://spaces.ac.cn/archives/7388
    """

    ds = D.shape[0]
    dsm = ds * ds
    A_eq = sp.dok_matrix((ds*2, dsm))
    
    for i in range(len(p)):
        
        A = sp.dok_matrix((ds, ds))
        A[i, :] = 1
        A_eq[i] = A.reshape(1, dsm) 
    
    for i in range(len(q)):
    
        A = sp.dok_matrix((ds, ds))
        A[:, i] = 1
        A_eq[i+ds] = A.reshape(1, dsm) 
        
    b_eq = np.concatenate([p, q])
    D = D.reshape(-1)
    result = linprog(D, A_eq=A_eq[:-1], b_eq=b_eq[:-1], method='interior-point', options={'sparse': True})
    
    return result.fun


def word_mover_distance(x, y):

    """
    Calculate word mover's distance (WMD)
    Adapted from: Jianlin Su @ https://spaces.ac.cn/archives/7388
    """

    p = np.ones(x.shape[0]) / x.shape[0]
    q = np.ones(y.shape[0]) / y.shape[0]
    D = np.sqrt(np.square(x[:, None] - y[None, :]).mean(axis=2))
    
    return wasserstein_distance(p, q, D)


def word_rotator_distance(x, y):

    """
    Calculate word rotator's distance (WRD)
    Adapted from: Jianlin Su @ https://spaces.ac.cn/archives/7388
    """
    
    x_norm = (x**2).sum(axis=1, keepdims=True)**0.5
    y_norm = (y**2).sum(axis=1, keepdims=True)**0.5
    p = x_norm[:, 0] / x_norm.sum()
    q = y_norm[:, 0] / y_norm.sum()
    D = 1 - np.dot(x / x_norm, (y / y_norm).T)
    
    return wasserstein_distance(p, q, D)


def spearmanr(x, y):

    """
    Compute Spearman rho rank-order correlation
    """

    x = np.column_stack((x, y))
    x_ranked = np.apply_along_axis(rankdata, 0, x)
    
    return np.corrcoef(x_ranked, rowvar=0)[0][1]


def permutation_corr(x, y, iter=1000):

    """
    Conduct permutation test for correlation coefficients
    Adapted from: Zitong Lu @ https://github.com/ZitongLu1996/NeuroRA/blob/master/neurora/stuff.py
    """

    rtest = spearmanr(x, y)

    ni = 0

    for i in range(iter):
            
        x_shuffle = np.random.permutation(x)
        y_shuffle = np.random.permutation(y)
        rperm = spearmanr(x_shuffle, y_shuffle)

        if rperm >= rtest:
            ni = ni + 1

    p = np.float64((ni+1)/(iter+1))

    return p


def load_mms_data():

    """
    Load surface-based multivariate morphometry statistics (MMS) data
    """

    order = [i + 102 for i in range(31) if i + 102 not in exclude_num]
    mms_la, mms_ra, mms_lh, mms_rh = [], [], [], []

    for n in order:

        f_la = open('./data/mms/amyg/lh/' + str(n) +
                    '_T1_LAmyg_60k_std_par_flowed_jfeat.m')
        f_ra = open('./data/mms/amyg/rh/' + str(n) +
                    '_T1_RAmyg_60k_std_par_flowed_jfeat.m')
        f_lh = open('./data/mms/hippo/lh/' + str(n) +
                    '_T1_LHippo_60k_std_par_flowed_jfeat.m')
        f_rh = open('./data/mms/hippo/rh/' + str(n) +
                    '_T1_RHippo_60k_std_par_flowed_jfeat.m')

        mms_la.append(
            [[float(j) for j in i.split('Jfeature=(')[-1][:-3].split(' ')[:4]]
             for i in f_la.readlines()[:15000]])
        mms_ra.append(
            [[float(j) for j in i.split('Jfeature=(')[-1][:-3].split(' ')[:4]]
             for i in f_ra.readlines()[:15000]])
        mms_lh.append(
            [[float(j) for j in i.split('Jfeature=(')[-1][:-3].split(' ')[:4]]
             for i in f_lh.readlines()[:15000]])
        mms_rh.append(
            [[float(j) for j in i.split('Jfeature=(')[-1][:-3].split(' ')[:4]]
             for i in f_rh.readlines()[:15000]])

        f_la.close()
        f_ra.close()
        f_lh.close()
        f_rh.close()

    return np.array(mms_la), np.array(mms_ra), np.array(mms_lh), np.array(
        mms_rh)


def load_rsfc_data():

    """
    Load Resting-state functional connectivity (rs-FC) data
    """

    f = open('./data/rsfc/rsfc.pkl', 'rb')
    rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh = pickle.load(f)

    return rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh


def load_IMQ_data(score='mvs', rank=True):

    """
    Load interactive mentalizing questionnaire (IMQ) data

    Parameters
    ----------
    score: string 'mvs' or 'uvs'. Default is 'mvs'.
        'Mvs' for multivariate representations, i.e., normalised item-wise IMQ scores, 
        'Uvs' for univariate representations, i.e, scalar summary IMQ scores.
    rank: bool True or False. Default is True.
        Assign ranks to data or not.

    See more info about IMQ: https://github.com/andlab-um/IMQ and https://doi.org/10.3389/fpsyg.2021.791835
    """

    if score == 'mvs':
        f = open('./data/imq/mvs.pkl', 'rb')
    if score == 'uvs':
        f = open('./data/imq/uvs.pkl', 'rb')

    SS, SO, OS = pickle.load(f)

    if score == 'uvs' and rank:
        return np.array(rankdata(SS)), np.array(rankdata(SO)), np.array(rankdata(OS))
    else:
        return np.array(SS), np.array(SO), np.array(OS)


def dissimilarity(x, y, cov=None, metric=None):

    """
    Calculate the dissimilarities

    Parameters
    ----------
    x: array.
        Vector x.
    y: array.
        Vector y.
    cov: array. Default is None.
        Covariance matrix for the computing of Mahalanobis distance.
    metric: string 'abs' or 'mean' or 'min' or 'abs*mean' or 'max' or 
            'pearson' or 'euclidean' or 'mahalanobis' or 'cosine' or 'manhattan' or
            'wmd' or 'wrd'. Default is None.
        The method to calculate the dissimilarities or distances.
        'Abs' for calculating absolute distance, 'mean' for calculating mean distance, 'min' for calculating minimum distance, 
        'abs*mean' for calculating the product of the absolute and minimum distance, 'max' for calculating maximum distance, 'pearson' for calculating Pearson distance, 
        'euclidean' for calculating Euclidean distance, 'mahalanobis' for calculating Mahalanobis distance, 'cosine' for calculating cosine distance, 
        'manhattan' for calculating Manhattan distance, 'wmd' for calculating word mover distance, 'wrd' for calculating word rotator distance. 
    """

    if metric in [
            'pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan'
    ]:
        x, y = x.reshape(-1), y.reshape(-1)

    if metric == 'abs':
        return abs(x-y)
    elif metric == 'mean':
        return 1-np.mean((x, y))/24
    elif metric == 'min':
        return 1-min(x, y)/24
    elif metric == 'abs*mean':
        return abs(x-y)*(1-np.mean((x, y))/24)
    elif metric == 'max':
        return max(x, y)/24
    elif metric == 'pearson':
        return distance.correlation(x, y)
    elif metric == 'euclidean':
        return np.linalg.norm(x - y)
    elif metric == 'mahalanobis':
        return np.sqrt(np.dot(np.dot(x-y, cov), x-y))
    elif metric == 'cosine':
        return 1-np.sum(x * y)/(np.linalg.norm(x) * np.linalg.norm(y))
    elif metric == 'manhattan':
        return np.sum(np.abs(x - y))
    elif metric == 'wmd':
        return word_mover_distance(x, y)
    elif metric == 'wrd':
        return word_rotator_distance(x, y)


def build_idm(data, excluded_pairs=None, metric=None):

    """
    Construct the inter-subject dissimilarity matrix (IDM)
    
    Parameters
    ----------
    data: array.
        It could be MMS or rs-FC or IMQ score.
    excluded_pairs: list. Default is None.
        Excluded non-informative subject pairs.
    metric: string 'abs' or 'mean' or 'min' or 'abs*mean' or 'max' or 
            'pearson' or 'euclidean' or 'mahalanobis' or 'cosine' or 'manhattan' or
            'wmd' or 'wrd'. Default is None.
        The method to calculate the dissimilarities or distances.
        'Abs' for calculating absolute distance, 'mean' for calculating mean distance, 'min' for calculating minimum distance, 
        'abs*mean' for calculating the product of the absolute and minimum distance, 'max' for calculating maximum distance, 'pearson' for calculating Pearson distance, 
        'euclidean' for calculating Euclidean distance, 'mahalanobis' for calculating Mahalanobis distance, 'cosine' for calculating cosine distance, 
        'manhattan' for calculating Manhattan distance, 'wmd' for calculating word mover distance, 'wrd' for calculating word rotator distance. 
    
    Returns
    -------
    result: array.
        The upper off-diagonal triangle of IDM (vector).
    """

    length = len(data)

    if metric == 'mahalanobis':
        cov = np.cov(data.reshape(length, -1), rowvar=False)
    else:
        cov = None

    if excluded_pairs != None:
        result = [
            dissimilarity(data[i], data[j], cov, metric)
            for i in range(length) for j in range(length) if i < j and (i, j) not in excluded_pairs 
        ]

    else:
        result = [
            dissimilarity(data[i], data[j], cov, metric)
            for i in range(length) for j in range(length) if i < j 
        ]

    result = np.array(result)

    return result


def patch(data, pm=None, ps=(10, 10), sd=None):

    """
    Patch the vertices and conduct global pooling operation within each patch for constructing MMS IDM
    
    Parameters
    ----------
    data: array.
        MMS data.
    pm: string 'max' or 'mean' or 'min' or 'max-mean' or 'max-min' or 'mean-min' or 'mmm'. Default is None.
        The method to conduct global pooling operation.
        'Max' for max-over-vertex pooling operation, 'mean' for mean-over-vertex pooling operation, 
        'min' for min-over-vertex pooling operation, 'max-mean' for the combination (i.e., concatenation) of max- and mean-over-vertex pooling operation,
        'max-min' for the combination (i.e., concatenation) of max- and min-over-vertex pooling operation,
        'mean-min' for the combination (i.e., concatenation) of mean- and min-over-vertex pooling operation,
        'mmm' for the combination (i.e., concatenation) of max-, mean- and min-over-vertex pooling operation.
    ps: tuple. Default is (10, 10).
        Patching size. it could be 5 × 5, 10 × 10, 25 × 25 and 50 × 50 vertices. 
    sd: string 'pearson' or 'euclidean' or 'mahalanobis' or 'cosine' or 'manhattan'. Default is None.
        The method to calculate the surface distance. 
        'pearson' for calculating Pearson distance, 'euclidean' for calculating Euclidean distance, 'mahalanobis' for calculating Mahalanobis distance, 
        'cosine' for calculating cosine distance, 'manhattan' for calculating Manhattan distance.

    Returns
    -------
    data_patch: array.
        The log-dimensional representation for high-dimensional MMS data.
    """

    data_patch = []
    for d in data:

        patch = []
        for i in range(0, 150, ps[1]):
            for j in range(0, 100, ps[0]):

                p = []
                for pi in range(ps[1]):
                    for pj in range(ps[0]):

                        p.append(d[(i + pi) * 100 + (j + pj)])

                if pm == None:
                    patch.append(np.array(p).reshape(-1))
                elif pm == 'max':
                    patch.append(np.max(np.array(p), axis=0))
                elif pm == 'mean':
                    patch.append(np.mean(np.array(p), axis=0))
                elif pm == 'min':
                    patch.append(np.min(np.array(p), axis=0))
                elif pm == 'max-mean':
                    patch.append(
                        np.array(
                            list(np.max(np.array(p), axis=0)) +
                            list(np.mean(np.array(p), axis=0))))
                elif pm == 'max-min':
                    patch.append(
                        np.array(
                            list(np.max(np.array(p), axis=0)) +
                            list(np.min(np.array(p), axis=0))))
                elif pm == 'mean-min':
                    patch.append(
                        np.array(
                            list(np.mean(np.array(p), axis=0)) +
                            list(np.min(np.array(p), axis=0))))
                elif pm == 'mmm':
                    patch.append(
                        np.array(
                            list(np.max(np.array(p), axis=0)) +
                            list(np.mean(np.array(p), axis=0)) +
                            list(np.min(np.array(p), axis=0))))

        if sd != None:
            patch = build_rdm(np.array(patch), metric=sd)

        data_patch.append(patch)

    data_patch = np.array(data_patch)

    return data_patch


def run(x, y, dsf_x, dsf_y, mms=None, iter=1000, idm=False):

    """
    Compare two IDMs using Spearman rho rank-order correlation

    Parameters
    ----------
    x: array.
        It could be rs-FC or IMQ score.
    y: array.
        It could be MMS or rs-FC.
    dsf_x: string 'abs' or 'mean' or 'min' or 'abs*mean' or 'max' or 
            'pearson' or 'euclidean' or 'mahalanobis' or 'cosine' or 'manhattan'.
        The method to calculate the dissimilarities or distances for the modality x.
        'Abs' for calculating absolute distance, 'mean' for calculating mean distance, 'min' for calculating minimum distance, 
        'abs*mean' for calculating the product of the absolute and minimum distance, 'max' for calculating maximum distance, 'pearson' for calculating Pearson distance, 
        'euclidean' for calculating Euclidean distance, 'mahalanobis' for calculating Mahalanobis distance, 'cosine' for calculating cosine distance, 
        'manhattan' for calculating Manhattan distance.
    dsf_y: string 'pearson' or 'euclidean' or 'mahalanobis' or 'cosine' or 'manhattan' or 'wmd' or 'wrd'.
        The method to calculate the dissimilarities or distances for the modality y.
        'Pearson' for calculating Pearson distance, 'euclidean' for calculating Euclidean distance, 'mahalanobis' for calculating Mahalanobis distance, 
        'cosine' for calculating cosine distance, 'manhattan' for calculating Manhattan distance, 'wmd' for calculating word mover distance, 'wrd' for calculating word rotator distance. 
    mms: tuple. Default is None.
        The parameter configuration for constructing MMS IDM.
        (P, PS, PM, SD), 'P': bool True or False. Conduct patching or not. 'PS', 'PM' and 'SD' refer to function 'patch'.
    iter: int. Default is 1000.
        The iteration numbers for the permutation test.
    idm: bool True or False. Default is False.
        If idm is True, then return Spearman rho rank-order correlation, x_idm and y_idm;
        if idm is False, then return Spearman rho rank-order correlation, p-value derived from the permutation test.
    """

    x_idm = build_idm(x, metric=dsf_x)

    if mms != None:

        P, PS, PM, SD = mms

        if P:
            y_idm = build_idm(patch(y, pm=PM, ps=PS, sd=SD), metric=dsf_y)

    else:
        y_idm = build_idm(y, metric=dsf_y)
        
    if idm:
        return spearmanr(x_idm, y_idm), x_idm, y_idm
    else:
        return spearmanr(x_idm, y_idm), permutation_corr(x_idm, y_idm, iter=iter)