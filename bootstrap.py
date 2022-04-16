# -*- coding: utf-8 -*-

__author__ = 'Zhaoning Li'

from util import *
import concurrent.futures


def IMQ_MMS_job(p):

    dsf_x, dsf_y, mms, data, SS, SO, OS, mms_la, mms_ra, mms_lh, mms_rh = p

    result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh = [], [], [], []
    result_so_mms_la, result_so_mms_ra, result_so_mms_lh, result_so_mms_rh = [], [], [], []
    result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh = [], [], [], []

    for i, order in enumerate(tqdm(data[0])):

        SS_sample = SS[order]
        SO_sample = SO[order]
        OS_sample = OS[order]

        mms_la_sample = mms_la[order]
        mms_ra_sample = mms_ra[order]
        mms_lh_sample = mms_lh[order]
        mms_rh_sample = mms_rh[order]

        SS_sample_idm = build_idm(SS_sample, excluded_pairs=data[1][i], metric=dsf_x)
        SO_sample_idm = build_idm(SO_sample, excluded_pairs=data[1][i], metric=dsf_x)
        OS_sample_idm = build_idm(OS_sample, excluded_pairs=data[1][i], metric=dsf_x)

        if mms != None:

            P, PS, PM, SD = mms
            mms_la_sample_idm = build_idm(patch(mms_la_sample, pm=PM, ps=PS, sd=SD), excluded_pairs=data[1][i], metric=dsf_y)
            mms_ra_sample_idm = build_idm(patch(mms_ra_sample, pm=PM, ps=PS, sd=SD), excluded_pairs=data[1][i], metric=dsf_y)
            mms_lh_sample_idm = build_idm(patch(mms_lh_sample, pm=PM, ps=PS, sd=SD), excluded_pairs=data[1][i], metric=dsf_y)
            mms_rh_sample_idm = build_idm(patch(mms_rh_sample, pm=PM, ps=PS, sd=SD), excluded_pairs=data[1][i], metric=dsf_y)

        else:

            P, PS, PM, SD = False, None, None, None
            mms_la_sample_idm = build_idm(mms_la_sample, excluded_pairs=data[1][i], metric=dsf_y)
            mms_ra_sample_idm = build_idm(mms_ra_sample, excluded_pairs=data[1][i], metric=dsf_y)
            mms_lh_sample_idm = build_idm(mms_lh_sample, excluded_pairs=data[1][i], metric=dsf_y)
            mms_rh_sample_idm = build_idm(mms_rh_sample, excluded_pairs=data[1][i], metric=dsf_y)

        result_ss_mms_la.append(
            spearmanr(SS_sample_idm, mms_la_sample_idm))
        result_ss_mms_ra.append(
            spearmanr(SS_sample_idm, mms_ra_sample_idm))
        result_ss_mms_lh.append(
            spearmanr(SS_sample_idm, mms_lh_sample_idm))
        result_ss_mms_rh.append(
            spearmanr(SS_sample_idm, mms_rh_sample_idm))

        result_so_mms_la.append(
            spearmanr(SO_sample_idm, mms_la_sample_idm))
        result_so_mms_ra.append(
            spearmanr(SO_sample_idm, mms_ra_sample_idm))
        result_so_mms_lh.append(
            spearmanr(SO_sample_idm, mms_lh_sample_idm))
        result_so_mms_rh.append(
            spearmanr(SO_sample_idm, mms_rh_sample_idm))

        result_os_mms_la.append(
            spearmanr(OS_sample_idm, mms_la_sample_idm))
        result_os_mms_ra.append(
            spearmanr(OS_sample_idm, mms_ra_sample_idm))
        result_os_mms_lh.append(
            spearmanr(OS_sample_idm, mms_lh_sample_idm))
        result_os_mms_rh.append(
            spearmanr(OS_sample_idm, mms_rh_sample_idm))

    configuration = dsf_x+'-'+dsf_y+'-' + \
        str(P)+'-'+str(PS)+'-'+str(PM)+'-'+str(SD)

    print(configuration)

    return configuration, result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra, result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh


def IMQ_rsFC_job(p):

    dsf_x, dsf_y, mms, data, SS, SO, OS, rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh = p

    result_ss_rsfc_la, result_ss_rsfc_ra, result_ss_rsfc_lh, result_ss_rsfc_rh = [], [], [], []
    result_so_rsfc_la, result_so_rsfc_ra, result_so_rsfc_lh, result_so_rsfc_rh = [], [], [], []
    result_os_rsfc_la, result_os_rsfc_ra, result_os_rsfc_lh, result_os_rsfc_rh = [], [], [], []

    for i, order in enumerate(tqdm(data[0])):

        SS_sample = SS[order]
        SO_sample = SO[order]
        OS_sample = OS[order]

        rsfc_la_sample = rsfc_la[order]
        rsfc_ra_sample = rsfc_ra[order]
        rsfc_lh_sample = rsfc_lh[order]
        rsfc_rh_sample = rsfc_rh[order]

        SS_sample_idm = build_idm(SS_sample, excluded_pairs=data[1][i], metric=dsf_x)
        SO_sample_idm = build_idm(SO_sample, excluded_pairs=data[1][i], metric=dsf_x)
        OS_sample_idm = build_idm(OS_sample, excluded_pairs=data[1][i], metric=dsf_x)

        rsfc_la_sample_idm = build_idm(rsfc_la_sample, excluded_pairs=data[1][i], metric=dsf_y)
        rsfc_ra_sample_idm = build_idm(rsfc_ra_sample, excluded_pairs=data[1][i], metric=dsf_y)
        rsfc_lh_sample_idm = build_idm(rsfc_lh_sample, excluded_pairs=data[1][i], metric=dsf_y)
        rsfc_rh_sample_idm = build_idm(rsfc_rh_sample, excluded_pairs=data[1][i], metric=dsf_y)

        result_ss_rsfc_la.append(
            spearmanr(SS_sample_idm, rsfc_la_sample_idm))
        result_ss_rsfc_ra.append(
            spearmanr(SS_sample_idm, rsfc_ra_sample_idm))
        result_ss_rsfc_lh.append(
            spearmanr(SS_sample_idm, rsfc_lh_sample_idm))
        result_ss_rsfc_rh.append(
            spearmanr(SS_sample_idm, rsfc_rh_sample_idm))

        result_so_rsfc_la.append(
            spearmanr(SO_sample_idm, rsfc_la_sample_idm))
        result_so_rsfc_ra.append(
            spearmanr(SO_sample_idm, rsfc_ra_sample_idm))
        result_so_rsfc_lh.append(
            spearmanr(SO_sample_idm, rsfc_lh_sample_idm))
        result_so_rsfc_rh.append(
            spearmanr(SO_sample_idm, rsfc_rh_sample_idm))

        result_os_rsfc_la.append(
            spearmanr(OS_sample_idm, rsfc_la_sample_idm))
        result_os_rsfc_ra.append(
            spearmanr(OS_sample_idm, rsfc_ra_sample_idm))
        result_os_rsfc_lh.append(
            spearmanr(OS_sample_idm, rsfc_lh_sample_idm))
        result_os_rsfc_rh.append(
            spearmanr(OS_sample_idm, rsfc_rh_sample_idm))

    configuration = dsf_x+'-'+dsf_y

    print(configuration)

    return configuration, result_ss_rsfc_la, result_ss_rsfc_ra, result_ss_rsfc_lh, result_ss_rsfc_rh, result_so_rsfc_la, result_so_rsfc_ra, result_so_rsfc_lh, result_so_rsfc_rh, result_os_rsfc_la, result_os_rsfc_ra, result_os_rsfc_lh, result_os_rsfc_rh


def rsFC_MMS_job(p):

    dsf_x, dsf_y, mms, data, rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh, mms_la, mms_ra, mms_lh, mms_rh = p

    result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh = [], [], [], []
    result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh = [], [], [], []
    result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh = [], [], [], []
    result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh = [], [], [], []

    for i, order in enumerate(tqdm(data[0])):

        rsfc_la_sample = rsfc_la[order]
        rsfc_ra_sample = rsfc_ra[order]
        rsfc_lh_sample = rsfc_lh[order]
        rsfc_rh_sample = rsfc_rh[order]

        mms_la_sample = mms_la[order]
        mms_ra_sample = mms_ra[order]
        mms_lh_sample = mms_lh[order]
        mms_rh_sample = mms_rh[order]

        rsfc_la_sample_idm = build_idm(rsfc_la_sample, excluded_pairs=data[1][i], metric=dsf_x)
        rsfc_ra_sample_idm = build_idm(rsfc_ra_sample, excluded_pairs=data[1][i], metric=dsf_x)
        rsfc_lh_sample_idm = build_idm(rsfc_lh_sample, excluded_pairs=data[1][i], metric=dsf_x)
        rsfc_rh_sample_idm = build_idm(rsfc_rh_sample, excluded_pairs=data[1][i], metric=dsf_x)

        if mms != None:

            P, PS, PM, SD = mms
            mms_la_sample_idm = build_idm(patch(mms_la_sample, pm=PM, ps=PS, sd=SD), excluded_pairs=data[1][i], metric=dsf_y)
            mms_ra_sample_idm = build_idm(patch(mms_ra_sample, pm=PM, ps=PS, sd=SD), excluded_pairs=data[1][i], metric=dsf_y)
            mms_lh_sample_idm = build_idm(patch(mms_lh_sample, pm=PM, ps=PS, sd=SD), excluded_pairs=data[1][i], metric=dsf_y)
            mms_rh_sample_idm = build_idm(patch(mms_rh_sample, pm=PM, ps=PS, sd=SD), excluded_pairs=data[1][i], metric=dsf_y)

        else:

            P, PS, PM, SD = False, None, None, None
            mms_la_sample_idm = build_idm(mms_la_sample, excluded_pairs=data[1][i], metric=dsf_y)
            mms_ra_sample_idm = build_idm(mms_ra_sample, excluded_pairs=data[1][i], metric=dsf_y)
            mms_lh_sample_idm = build_idm(mms_lh_sample, excluded_pairs=data[1][i], metric=dsf_y)
            mms_rh_sample_idm = build_idm(mms_rh_sample, excluded_pairs=data[1][i], metric=dsf_y)

        result_rsfc_la_mms_la.append(
            spearmanr(rsfc_la_sample_idm, mms_la_sample_idm))
        result_rsfc_la_mms_ra.append(
            spearmanr(rsfc_la_sample_idm, mms_ra_sample_idm))
        result_rsfc_la_mms_lh.append(
            spearmanr(rsfc_la_sample_idm, mms_lh_sample_idm))
        result_rsfc_la_mms_rh.append(
            spearmanr(rsfc_la_sample_idm, mms_rh_sample_idm))

        result_rsfc_ra_mms_la.append(
            spearmanr(rsfc_ra_sample_idm, mms_la_sample_idm))
        result_rsfc_ra_mms_ra.append(
            spearmanr(rsfc_ra_sample_idm, mms_ra_sample_idm))
        result_rsfc_ra_mms_lh.append(
            spearmanr(rsfc_ra_sample_idm, mms_lh_sample_idm))
        result_rsfc_ra_mms_rh.append(
            spearmanr(rsfc_ra_sample_idm, mms_rh_sample_idm))

        result_rsfc_lh_mms_la.append(
            spearmanr(rsfc_lh_sample_idm, mms_la_sample_idm))
        result_rsfc_lh_mms_ra.append(
            spearmanr(rsfc_lh_sample_idm, mms_ra_sample_idm))
        result_rsfc_lh_mms_lh.append(
            spearmanr(rsfc_lh_sample_idm, mms_lh_sample_idm))
        result_rsfc_lh_mms_rh.append(
            spearmanr(rsfc_lh_sample_idm, mms_rh_sample_idm))

        result_rsfc_rh_mms_la.append(
            spearmanr(rsfc_rh_sample_idm, mms_la_sample_idm))
        result_rsfc_rh_mms_ra.append(
            spearmanr(rsfc_rh_sample_idm, mms_ra_sample_idm))
        result_rsfc_rh_mms_lh.append(
            spearmanr(rsfc_rh_sample_idm, mms_lh_sample_idm))
        result_rsfc_rh_mms_rh.append(
            spearmanr(rsfc_rh_sample_idm, mms_rh_sample_idm))

    configuration = dsf_x+'-'+dsf_y+'-' + \
        str(P)+'-'+str(PS)+'-'+str(PM)+'-'+str(SD)

    print(configuration)

    return configuration, result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh, result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh, result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh, result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh


def bootstrapping(mode, score='mvs', mms_pc=None, mw=48, size=1000):

    """
    Conduct subject-wise bootstapping for possible parameter configurations

    Parameters
    ----------
    mode: string 'IMQ-MMS' or 'IMQ-rsFC' or 'rsFC-MMS'.
        The combination between two modalities.
    score: string 'mvs' or 'uvs'. Default is 'mvs'.
        'Mvs' for multivariate representations, i.e., normalised item-wise IMQ scores, 
        'Uvs' for univariate representations, i.e, scalar summary IMQ scores.
    mms_pc: string 'NPP' or 'PP' or 'PPw' or 'CPP-SD'. Default is None.
        The pipeline to construct inter-subject MMS dissimilarity matrix.
        'NPP' for without patching and pooling operations and computing surface distance,
        'PP' or 'PPw' (wmd and wrd) for with patching and pooling operations but without computing surface distance,
        'CPP-SD' for computing patching and pooling operations-based surface distance.
    mw: int. Default is 48.
        Maximum number of workers that can be simulatenously active at the same time
    size: int. Default is 1000.
        The repetitions numbers for subject-wise bootstapping.
    """

    parameter = []
    order = np.arange(24)
    data = []

    for i in range(size):

        np.random.seed(i)
        data.append(np.random.choice(order, size=24))

    del_data = [excluded_pair(d) for d in data]
    data = [data, del_data]

    if mode == 'IMQ-MMS':

        SS, SO, OS = load_IMQ_data(score=score)
        mms_la, mms_ra, mms_lh, mms_rh = load_mms_data()
        result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh = [], [], [], []
        result_so_mms_la, result_so_mms_ra, result_so_mms_lh, result_so_mms_rh = [], [], [], []
        result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh = [], [], [], []

        if score == 'mvs':
            DSF_X = ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']
        if score == 'uvs':
            DSF_X = ['abs', 'mean', 'min', 'abs*mean', 'max']

        if mms_pc == 'NPP':

            for dsf_x in DSF_X:
                for dsf_y in ['pearson', 'euclidean', 'cosine', 'manhattan']:

                    mms = None
                    parameter.append(
                        (dsf_x, dsf_y, mms, data, SS, SO, OS, mms_la, mms_ra, mms_lh, mms_rh))

            with concurrent.futures.ProcessPoolExecutor(max_workers=mw) as executor:
                results = list(
                    tqdm(executor.map(IMQ_MMS_job, parameter), total=5*4))

        elif mms_pc == 'PP':

            for dsf_x in DSF_X:
                for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                    for P in [True]:
                        for PS in [(5, 5), (10, 10), (25, 25), (50, 50)]:
                            for PM in ['max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                                for SD in [None]:

                                    mms = P, PS, PM, SD
                                    parameter.append(
                                        (dsf_x, dsf_y, mms, data, SS, SO, OS, mms_la, mms_ra, mms_lh, mms_rh))

            with concurrent.futures.ProcessPoolExecutor(max_workers=mw) as executor:
                results = list(
                    tqdm(executor.map(IMQ_MMS_job, parameter), total=5*5*4*7))

        elif mms_pc == 'PPw':

            for dsf_x in DSF_X:
                for dsf_y in ['wrd', 'wmd']:
                    for P in [True]:
                        for PS in [(25, 25), (50, 50)]:
                            for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                                for SD in [None]:

                                    mms = P, PS, PM, SD
                                    parameter.append(
                                        (dsf_x, dsf_y, mms, data, SS, SO, OS, mms_la, mms_ra, mms_lh, mms_rh))

            with concurrent.futures.ProcessPoolExecutor(max_workers=mw) as executor:
                results = list(
                    tqdm(executor.map(IMQ_MMS_job, parameter), total=5*2*2*8))

        elif mms_pc == 'CPP-SD':

            for dsf_x in DSF_X:
                for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                    for P in [True]:
                        for PS in [(25, 25), (50, 50)]:
                            for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                                for SD in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:

                                    mms = P, PS, PM, SD
                                    parameter.append(
                                        (dsf_x, dsf_y, mms, data, SS, SO, OS, mms_la, mms_ra, mms_lh, mms_rh))

            with concurrent.futures.ProcessPoolExecutor(max_workers=mw) as executor:
                results = list(
                    tqdm(executor.map(IMQ_MMS_job, parameter), total=5*5*2*8*5))

        result_ss_mms_la = [(r[0], r[1]) for r in results]
        result_ss_mms_ra = [(r[0], r[2]) for r in results]
        result_ss_mms_lh = [(r[0], r[3]) for r in results]
        result_ss_mms_rh = [(r[0], r[4]) for r in results]

        result_so_mms_la = [(r[0], r[5]) for r in results]
        result_so_mms_ra = [(r[0], r[6]) for r in results]
        result_so_mms_lh = [(r[0], r[7]) for r in results]
        result_so_mms_rh = [(r[0], r[8]) for r in results]

        result_os_mms_la = [(r[0], r[9]) for r in results]
        result_os_mms_ra = [(r[0], r[10]) for r in results]
        result_os_mms_lh = [(r[0], r[11]) for r in results]
        result_os_mms_rh = [(r[0], r[12]) for r in results]

        f = open('./bs_log/IMQ-MMS/'+score+'/'+mms_pc+'/results.pkl', 'wb')
        pickle.dump((result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra,
                    result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh), f)

    elif mode == 'IMQ-rsFC':

        SS, SO, OS = load_IMQ_data(score=score)
        rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh = load_rsfc_data()
        result_ss_rsfc_la, result_ss_rsfc_ra, result_ss_rsfc_lh, result_ss_rsfc_rh = [], [], [], []
        result_so_rsfc_la, result_so_rsfc_ra, result_so_rsfc_lh, result_so_rsfc_rh = [], [], [], []
        result_os_rsfc_la, result_os_rsfc_ra, result_os_rsfc_lh, result_os_rsfc_rh = [], [], [], []

        if score == 'mvs':
            DSF_X = ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']
        if score == 'uvs':
            DSF_X = ['abs', 'mean', 'min', 'abs*mean', 'max']

        for dsf_x in DSF_X:
            for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:

                mms = None
                parameter.append((dsf_x, dsf_y, mms, data, SS, SO, OS, rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh))

        with concurrent.futures.ProcessPoolExecutor(max_workers=mw) as executor:
            results = list(
                tqdm(executor.map(IMQ_rsFC_job, parameter), total=5*5))

        result_ss_rsfc_la = [(r[0], r[1]) for r in results]
        result_ss_rsfc_ra = [(r[0], r[2]) for r in results]
        result_ss_rsfc_lh = [(r[0], r[3]) for r in results]
        result_ss_rsfc_rh = [(r[0], r[4]) for r in results]

        result_so_rsfc_la = [(r[0], r[5]) for r in results]
        result_so_rsfc_ra = [(r[0], r[6]) for r in results]
        result_so_rsfc_lh = [(r[0], r[7]) for r in results]
        result_so_rsfc_rh = [(r[0], r[8]) for r in results]

        result_os_rsfc_la = [(r[0], r[9]) for r in results]
        result_os_rsfc_ra = [(r[0], r[10]) for r in results]
        result_os_rsfc_lh = [(r[0], r[11]) for r in results]
        result_os_rsfc_rh = [(r[0], r[12]) for r in results]

        f = open('./bs_log/IMQ-rsFC/'+score+'/results.pkl', 'wb')
        pickle.dump((result_ss_rsfc_la, result_ss_rsfc_ra, result_ss_rsfc_lh, result_ss_rsfc_rh, result_so_rsfc_la, result_so_rsfc_ra,
                    result_so_rsfc_lh, result_so_rsfc_rh, result_os_rsfc_la, result_os_rsfc_ra, result_os_rsfc_lh, result_os_rsfc_rh), f)

    elif mode == 'rsFC-MMS':

        rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh = load_rsfc_data()
        mms_la, mms_ra, mms_lh, mms_rh = load_mms_data()
        result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh = [], [], [], []
        result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh = [], [], [], []
        result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh = [], [], [], []
        result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh = [], [], [], []

        if mms_pc == 'NPP':

            for dsf_x in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for dsf_y in ['pearson', 'euclidean', 'cosine', 'manhattan']:

                    mms = None
                    parameter.append((dsf_x, dsf_y, mms, data, rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh, mms_la, mms_ra, mms_lh, mms_rh))

            with concurrent.futures.ProcessPoolExecutor(max_workers=mw) as executor:
                results = list(
                    tqdm(executor.map(rsFC_MMS_job, parameter), total=5*4))

        elif mms_pc == 'PP':

            for dsf_x in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                    for P in [True]:
                        for PS in [(5, 5), (10, 10), (25, 25), (50, 50)]:
                            for PM in ['max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                                for SD in [None]:

                                    mms = P, PS, PM, SD
                                    parameter.append((dsf_x, dsf_y, mms, data, rsfc_la,
                                                     rsfc_ra, rsfc_lh, rsfc_rh, mms_la, mms_ra, mms_lh, mms_rh))

            with concurrent.futures.ProcessPoolExecutor(max_workers=mw) as executor:
                results = list(
                    tqdm(executor.map(rsFC_MMS_job, parameter), total=5*5*4*7))

        elif mms_pc == 'PPw':

            for dsf_x in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for dsf_y in ['wrd', 'wmd']:
                    for P in [True]:
                        for PS in [(25, 25), (50, 50)]:
                            for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                                for SD in [None]:

                                    mms = P, PS, PM, SD
                                    parameter.append(
                                        (dsf_x, dsf_y, mms, data, rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh, mms_la, mms_ra, mms_lh, mms_rh))

            with concurrent.futures.ProcessPoolExecutor(max_workers=mw) as executor:
                results = list(
                    tqdm(executor.map(rsFC_MMS_job, parameter), total=5*2*2*8))

        elif mms_pc == 'CPP-SD':

            for dsf_x in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                    for P in [True]:
                        for PS in [(25, 25), (50, 50)]:
                            for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                                for SD in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:

                                    mms = P, PS, PM, SD
                                    parameter.append((dsf_x, dsf_y, mms, data, rsfc_la,
                                                     rsfc_ra, rsfc_lh, rsfc_rh, mms_la, mms_ra, mms_lh, mms_rh))

            with concurrent.futures.ProcessPoolExecutor(max_workers=mw) as executor:
                results = list(
                    tqdm(executor.map(rsFC_MMS_job, parameter), total=5*5*2*8*5))

        result_rsfc_la_mms_la = [(r[0], r[1]) for r in results]
        result_rsfc_la_mms_ra = [(r[0], r[2]) for r in results]
        result_rsfc_la_mms_lh = [(r[0], r[3]) for r in results]
        result_rsfc_la_mms_rh = [(r[0], r[4]) for r in results]

        result_rsfc_ra_mms_la = [(r[0], r[5]) for r in results]
        result_rsfc_ra_mms_ra = [(r[0], r[6]) for r in results]
        result_rsfc_ra_mms_lh = [(r[0], r[7]) for r in results]
        result_rsfc_ra_mms_rh = [(r[0], r[8]) for r in results]

        result_rsfc_lh_mms_la = [(r[0], r[9]) for r in results]
        result_rsfc_lh_mms_ra = [(r[0], r[10]) for r in results]
        result_rsfc_lh_mms_lh = [(r[0], r[11]) for r in results]
        result_rsfc_lh_mms_rh = [(r[0], r[12]) for r in results]

        result_rsfc_rh_mms_la = [(r[0], r[13]) for r in results]
        result_rsfc_rh_mms_ra = [(r[0], r[14]) for r in results]
        result_rsfc_rh_mms_lh = [(r[0], r[15]) for r in results]
        result_rsfc_rh_mms_rh = [(r[0], r[16]) for r in results]

        f = open('./bs_log/rsFC-MMS/'+mms_pc+'/results.pkl', 'wb')
        pickle.dump((result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh, result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh,
                    result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh, result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh), f)


if __name__ == '__main__':

    #mw = 10
    #size = 1000
    #score = 'mvs'
    
    #bootstrapping('IMQ-MMS', score, mms_pc='NPP', mw=mw, size=size)
    #bootstrapping('IMQ-MMS', score, mms_pc='PP', mw=mw, size=size)
    #bootstrapping('IMQ-MMS', score, mms_pc='PPw', mw=mw, size=size)
    #bootstrapping('IMQ-MMS', score, mms_pc='CPP-SD', mw=mw, size=size)
    #bootstrapping('IMQ-rsFC', score, mw=mw, size=size)
    #bootstrapping('rsFC-MMS', mms_pc='NPP', mw=mw, size=size)
    #bootstrapping('rsFC-MMS', mms_pc='PP', mw=mw, size=size)
    #bootstrapping('rsFC-MMS', mms_pc='PPw', mw=mw, size=size)
    #bootstrapping('rsFC-MMS', mms_pc='CPP-SD', mw=mw, size=size)

    #score = 'uvs'
    
    #bootstrapping('IMQ-MMS', score, mms_pc='NPP', mw=mw, size=size)
    #bootstrapping('IMQ-MMS', score, mms_pc='PP', mw=mw, size=size)
    #bootstrapping('IMQ-MMS', score, mms_pc='PPw', mw=mw, size=size)
    #bootstrapping('IMQ-MMS', score, mms_pc='CPP-SD', mw=mw, size=size)
    #bootstrapping('IMQ-rsFC', score, mw=mw, size=size)