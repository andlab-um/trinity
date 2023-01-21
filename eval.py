# -*- coding: utf-8 -*-

__author__ = 'Zhaoning Li'

from util import *
import scipy.stats as st
from statsmodels.stats.multitest import multipletests
from scipy.stats import kstest, wilcoxon
from pandas.core.frame import DataFrame
from scipy.stats import zscore
from pymer4.models import Lmer
import pandas as pd

SS_mvs, SO_mvs, OS_mvs = load_IMQ_data(score='mvs')
SS_uvs, SO_uvs, OS_uvs = load_IMQ_data(score='uvs')
mms_la, mms_ra, mms_lh, mms_rh = load_mms_data()
rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh = load_rsfc_data()

DSF_X_mvs = ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']
DSF_X_uvs = ['abs', 'mean', 'min', 'abs*mean', 'max']


def load_bootstrapping_results(mode):
    
    """
    Load bootstrapping reuslts from the folder 'bs_log'

    Parameters
    ----------
    mode: string 'IMQ-MMS' or 'IMQ-rsFC' or 'rsFC-MMS'.
        The combination between two modalities.
    """

    if mode == 'IMQ-MMS':

        f_np_mvs = open('./bs_log/IMQ-MMS/mvs/NPP/results.pkl', 'rb')
        result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra, result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh = pickle.load(
            f_np_mvs)

        results_np_mvs = (result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra,
                          result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh)

        f_np_uvs = open('./bs_log/IMQ-MMS/uvs/NPP/results.pkl', 'rb')
        result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra, result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh = pickle.load(
            f_np_uvs)

        results_np_uvs = (result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra,
                          result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh)

        results_npp = [results_np_mvs[i]+results_np_uvs[i] for i in range(12)]

        f_p_mvs = open('./bs_log/IMQ-MMS/mvs/PP/results.pkl', 'rb')
        result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra, result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh = pickle.load(
            f_p_mvs)

        results_p_mvs = (result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra,
                         result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh)

        f_p_uvs = open('./bs_log/IMQ-MMS/uvs/PP/results.pkl', 'rb')
        result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra, result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh = pickle.load(
            f_p_uvs)

        results_p_uvs = (result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra,
                         result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh)

        results_pp = [results_p_mvs[i]+results_p_uvs[i] for i in range(12)]

        f_pw_mvs = open('./bs_log/IMQ-MMS/mvs/PPw/results.pkl', 'rb')
        result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra, result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh = pickle.load(
            f_pw_mvs)

        results_pw_mvs = (result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra,
                          result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh)

        f_pw_uvs = open('./bs_log/IMQ-MMS/uvs/PPw/results.pkl', 'rb')
        result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra, result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh = pickle.load(
            f_pw_uvs)

        results_pw_uvs = (result_ss_mms_la, result_ss_mms_ra, result_ss_mms_lh, result_ss_mms_rh, result_so_mms_la, result_so_mms_ra,
                          result_so_mms_lh, result_so_mms_rh, result_os_mms_la, result_os_mms_ra, result_os_mms_lh, result_os_mms_rh)

        results_ppw = [results_pw_mvs[i]+results_pw_uvs[i] for i in range(12)]

        f_psc_mvs_1 = open('./bs_log/IMQ-MMS/mvs/CPP-SD/results_1.pkl', 'rb')
        result_ss_mms_la_mvs_1, result_ss_mms_ra_mvs_1, result_ss_mms_lh_mvs_1, result_ss_mms_rh_mvs_1, \
        result_so_mms_la_mvs_1, result_so_mms_ra_mvs_1, result_so_mms_lh_mvs_1, result_so_mms_rh_mvs_1, \
        result_os_mms_la_mvs_1, result_os_mms_ra_mvs_1, result_os_mms_lh_mvs_1, result_os_mms_rh_mvs_1 = pickle.load(
                f_psc_mvs_1)

        f_psc_mvs_2 = open('./bs_log/IMQ-MMS/mvs/CPP-SD/results_2.pkl', 'rb')
        result_ss_mms_la_mvs_2, result_ss_mms_ra_mvs_2, result_ss_mms_lh_mvs_2, result_ss_mms_rh_mvs_2, \
        result_so_mms_la_mvs_2, result_so_mms_ra_mvs_2, result_so_mms_lh_mvs_2, result_so_mms_rh_mvs_2, \
        result_os_mms_la_mvs_2, result_os_mms_ra_mvs_2, result_os_mms_lh_mvs_2, result_os_mms_rh_mvs_2 = pickle.load(
                f_psc_mvs_2)

        f_psc_mvs_3 = open('./bs_log/IMQ-MMS/mvs/CPP-SD/results_3.pkl', 'rb')
        result_ss_mms_la_mvs_3, result_ss_mms_ra_mvs_3, result_ss_mms_lh_mvs_3, result_ss_mms_rh_mvs_3, \
        result_so_mms_la_mvs_3, result_so_mms_ra_mvs_3, result_so_mms_lh_mvs_3, result_so_mms_rh_mvs_3, \
        result_os_mms_la_mvs_3, result_os_mms_ra_mvs_3, result_os_mms_lh_mvs_3, result_os_mms_rh_mvs_3 = pickle.load(
                f_psc_mvs_3)

        f_psc_mvs_4_1 = open('./bs_log/IMQ-MMS/mvs/CPP-SD/results_4_1.pkl', 'rb')
        result_ss_mms_la_mvs_4_1, result_ss_mms_ra_mvs_4_1, result_ss_mms_lh_mvs_4_1, result_ss_mms_rh_mvs_4_1, \
        result_so_mms_la_mvs_4_1, result_so_mms_ra_mvs_4_1, result_so_mms_lh_mvs_4_1, result_so_mms_rh_mvs_4_1, \
        result_os_mms_la_mvs_4_1, result_os_mms_ra_mvs_4_1, result_os_mms_lh_mvs_4_1, result_os_mms_rh_mvs_4_1 = pickle.load(
                f_psc_mvs_4_1)

        f_psc_mvs_4_2 = open('./bs_log/IMQ-MMS/mvs/CPP-SD/results_4_2.pkl', 'rb')
        result_ss_mms_la_mvs_4_2, result_ss_mms_ra_mvs_4_2, result_ss_mms_lh_mvs_4_2, result_ss_mms_rh_mvs_4_2, \
        result_so_mms_la_mvs_4_2, result_so_mms_ra_mvs_4_2, result_so_mms_lh_mvs_4_2, result_so_mms_rh_mvs_4_2, \
        result_os_mms_la_mvs_4_2, result_os_mms_ra_mvs_4_2, result_os_mms_lh_mvs_4_2, result_os_mms_rh_mvs_4_2 = pickle.load(
                f_psc_mvs_4_2)

        f_psc_mvs_5 = open('./bs_log/IMQ-MMS/mvs/CPP-SD/results_5.pkl', 'rb')
        result_ss_mms_la_mvs_5, result_ss_mms_ra_mvs_5, result_ss_mms_lh_mvs_5, result_ss_mms_rh_mvs_5, \
        result_so_mms_la_mvs_5, result_so_mms_ra_mvs_5, result_so_mms_lh_mvs_5, result_so_mms_rh_mvs_5, \
        result_os_mms_la_mvs_5, result_os_mms_ra_mvs_5, result_os_mms_lh_mvs_5, result_os_mms_rh_mvs_5 = pickle.load(
                f_psc_mvs_5)

        results_psc_mvs = (
            result_ss_mms_la_mvs_1+result_ss_mms_la_mvs_2+result_ss_mms_la_mvs_3+result_ss_mms_la_mvs_4_1+result_ss_mms_la_mvs_4_2+result_ss_mms_la_mvs_5,
            result_ss_mms_ra_mvs_1+result_ss_mms_ra_mvs_2+result_ss_mms_ra_mvs_3+result_ss_mms_ra_mvs_4_1+result_ss_mms_ra_mvs_4_2+result_ss_mms_ra_mvs_5,
            result_ss_mms_lh_mvs_1+result_ss_mms_lh_mvs_2+result_ss_mms_lh_mvs_3+result_ss_mms_lh_mvs_4_1+result_ss_mms_lh_mvs_4_2+result_ss_mms_lh_mvs_5,
            result_ss_mms_rh_mvs_1+result_ss_mms_rh_mvs_2+result_ss_mms_rh_mvs_3+result_ss_mms_rh_mvs_4_1+result_ss_mms_rh_mvs_4_2+result_ss_mms_rh_mvs_5,
            result_so_mms_la_mvs_1+result_so_mms_la_mvs_2+result_so_mms_la_mvs_3+result_so_mms_la_mvs_4_1+result_so_mms_la_mvs_4_2+result_so_mms_la_mvs_5,
            result_so_mms_ra_mvs_1+result_so_mms_ra_mvs_2+result_so_mms_ra_mvs_3+result_so_mms_ra_mvs_4_1+result_so_mms_ra_mvs_4_2+result_so_mms_ra_mvs_5,
            result_so_mms_lh_mvs_1+result_so_mms_lh_mvs_2+result_so_mms_lh_mvs_3+result_so_mms_lh_mvs_4_1+result_so_mms_lh_mvs_4_2+result_so_mms_lh_mvs_5,
            result_so_mms_rh_mvs_1+result_so_mms_rh_mvs_2+result_so_mms_rh_mvs_3+result_so_mms_rh_mvs_4_1+result_so_mms_rh_mvs_4_2+result_so_mms_rh_mvs_5,
            result_os_mms_la_mvs_1+result_os_mms_la_mvs_2+result_os_mms_la_mvs_3+result_os_mms_la_mvs_4_1+result_os_mms_la_mvs_4_2+result_os_mms_la_mvs_5,
            result_os_mms_ra_mvs_1+result_os_mms_ra_mvs_2+result_os_mms_ra_mvs_3+result_os_mms_ra_mvs_4_1+result_os_mms_ra_mvs_4_2+result_os_mms_ra_mvs_5,
            result_os_mms_lh_mvs_1+result_os_mms_lh_mvs_2+result_os_mms_lh_mvs_3+result_os_mms_lh_mvs_4_1+result_os_mms_lh_mvs_4_2+result_os_mms_lh_mvs_5,
            result_os_mms_rh_mvs_1+result_os_mms_rh_mvs_2+result_os_mms_rh_mvs_3+result_os_mms_rh_mvs_4_1+result_os_mms_rh_mvs_4_2+result_os_mms_rh_mvs_5)

        f_psc_uvs_1 = open('./bs_log/IMQ-MMS/uvs/CPP-SD/results_1.pkl', 'rb')
        result_ss_mms_la_uvs_1, result_ss_mms_ra_uvs_1, result_ss_mms_lh_uvs_1, result_ss_mms_rh_uvs_1, \
            result_so_mms_la_uvs_1, result_so_mms_ra_uvs_1, result_so_mms_lh_uvs_1, result_so_mms_rh_uvs_1, \
            result_os_mms_la_uvs_1, result_os_mms_ra_uvs_1, result_os_mms_lh_uvs_1, result_os_mms_rh_uvs_1 = pickle.load(
                f_psc_uvs_1)

        f_psc_uvs_2 = open('./bs_log/IMQ-MMS/uvs/CPP-SD/results_2.pkl', 'rb')
        result_ss_mms_la_uvs_2, result_ss_mms_ra_uvs_2, result_ss_mms_lh_uvs_2, result_ss_mms_rh_uvs_2, \
            result_so_mms_la_uvs_2, result_so_mms_ra_uvs_2, result_so_mms_lh_uvs_2, result_so_mms_rh_uvs_2, \
            result_os_mms_la_uvs_2, result_os_mms_ra_uvs_2, result_os_mms_lh_uvs_2, result_os_mms_rh_uvs_2 = pickle.load(
                f_psc_uvs_2)

        f_psc_uvs_3 = open('./bs_log/IMQ-MMS/uvs/CPP-SD/results_3.pkl', 'rb')
        result_ss_mms_la_uvs_3, result_ss_mms_ra_uvs_3, result_ss_mms_lh_uvs_3, result_ss_mms_rh_uvs_3, \
            result_so_mms_la_uvs_3, result_so_mms_ra_uvs_3, result_so_mms_lh_uvs_3, result_so_mms_rh_uvs_3, \
            result_os_mms_la_uvs_3, result_os_mms_ra_uvs_3, result_os_mms_lh_uvs_3, result_os_mms_rh_uvs_3 = pickle.load(
                f_psc_uvs_3)

        f_psc_uvs_4 = open('./bs_log/IMQ-MMS/uvs/CPP-SD/results_4.pkl', 'rb')
        result_ss_mms_la_uvs_4, result_ss_mms_ra_uvs_4, result_ss_mms_lh_uvs_4, result_ss_mms_rh_uvs_4, \
            result_so_mms_la_uvs_4, result_so_mms_ra_uvs_4, result_so_mms_lh_uvs_4, result_so_mms_rh_uvs_4, \
            result_os_mms_la_uvs_4, result_os_mms_ra_uvs_4, result_os_mms_lh_uvs_4, result_os_mms_rh_uvs_4 = pickle.load(
                f_psc_uvs_4)

        f_psc_uvs_5 = open('./bs_log/IMQ-MMS/uvs/CPP-SD/results_5.pkl', 'rb')
        result_ss_mms_la_uvs_5, result_ss_mms_ra_uvs_5, result_ss_mms_lh_uvs_5, result_ss_mms_rh_uvs_5, \
            result_so_mms_la_uvs_5, result_so_mms_ra_uvs_5, result_so_mms_lh_uvs_5, result_so_mms_rh_uvs_5, \
            result_os_mms_la_uvs_5, result_os_mms_ra_uvs_5, result_os_mms_lh_uvs_5, result_os_mms_rh_uvs_5 = pickle.load(
                f_psc_uvs_5)

        results_psc_uvs = (
            result_ss_mms_la_uvs_1+result_ss_mms_la_uvs_2+result_ss_mms_la_uvs_3+result_ss_mms_la_uvs_4+result_ss_mms_la_uvs_5,
            result_ss_mms_ra_uvs_1+result_ss_mms_ra_uvs_2+result_ss_mms_ra_uvs_3+result_ss_mms_ra_uvs_4+result_ss_mms_ra_uvs_5,
            result_ss_mms_lh_uvs_1+result_ss_mms_lh_uvs_2+result_ss_mms_lh_uvs_3+result_ss_mms_lh_uvs_4+result_ss_mms_lh_uvs_5,
            result_ss_mms_rh_uvs_1+result_ss_mms_rh_uvs_2+result_ss_mms_rh_uvs_3+result_ss_mms_rh_uvs_4+result_ss_mms_rh_uvs_5,
            result_so_mms_la_uvs_1+result_so_mms_la_uvs_2+result_so_mms_la_uvs_3+result_so_mms_la_uvs_4+result_so_mms_la_uvs_5,
            result_so_mms_ra_uvs_1+result_so_mms_ra_uvs_2+result_so_mms_ra_uvs_3+result_so_mms_ra_uvs_4+result_so_mms_ra_uvs_5,
            result_so_mms_lh_uvs_1+result_so_mms_lh_uvs_2+result_so_mms_lh_uvs_3+result_so_mms_lh_uvs_4+result_so_mms_lh_uvs_5,
            result_so_mms_rh_uvs_1+result_so_mms_rh_uvs_2+result_so_mms_rh_uvs_3+result_so_mms_rh_uvs_4+result_so_mms_rh_uvs_5,
            result_os_mms_la_uvs_1+result_os_mms_la_uvs_2+result_os_mms_la_uvs_3+result_os_mms_la_uvs_4+result_os_mms_la_uvs_5,
            result_os_mms_ra_uvs_1+result_os_mms_ra_uvs_2+result_os_mms_ra_uvs_3+result_os_mms_ra_uvs_4+result_os_mms_ra_uvs_5,
            result_os_mms_lh_uvs_1+result_os_mms_lh_uvs_2+result_os_mms_lh_uvs_3+result_os_mms_lh_uvs_4+result_os_mms_lh_uvs_5,
            result_os_mms_rh_uvs_1+result_os_mms_rh_uvs_2+result_os_mms_rh_uvs_3+result_os_mms_rh_uvs_4+result_os_mms_rh_uvs_5
        )

        results_cpp_sd = [results_psc_mvs[i]+results_psc_uvs[i] for i in range(12)]

        return results_npp, results_pp, results_ppw, results_cpp_sd

    elif mode == 'IMQ-rsFC':

        f_mvs = open('./bs_log/IMQ-rsFC/mvs/results.pkl', 'rb')
        result_ss_rsfc_la_mvs, result_ss_rsfc_ra_mvs, result_ss_rsfc_lh_mvs, result_ss_rsfc_rh_mvs, \
        result_so_rsfc_la_mvs, result_so_rsfc_ra_mvs, result_so_rsfc_lh_mvs, result_so_rsfc_rh_mvs, \
        result_os_rsfc_la_mvs, result_os_rsfc_ra_mvs, result_os_rsfc_lh_mvs, result_os_rsfc_rh_mvs = pickle.load(
            f_mvs)

        f_uvs = open('./bs_log/IMQ-rsFC/uvs/results.pkl', 'rb')
        result_ss_rsfc_la_uvs, result_ss_rsfc_ra_uvs, result_ss_rsfc_lh_uvs, result_ss_rsfc_rh_uvs, \
        result_so_rsfc_la_uvs, result_so_rsfc_ra_uvs, result_so_rsfc_lh_uvs, result_so_rsfc_rh_uvs, \
        result_os_rsfc_la_uvs, result_os_rsfc_ra_uvs, result_os_rsfc_lh_uvs, result_os_rsfc_rh_uvs = pickle.load(
        f_uvs)

        results = (
            result_ss_rsfc_la_mvs+result_ss_rsfc_la_uvs,
            result_ss_rsfc_ra_mvs+result_ss_rsfc_ra_uvs,
            result_ss_rsfc_lh_mvs+result_ss_rsfc_lh_uvs,
            result_ss_rsfc_rh_mvs+result_ss_rsfc_rh_uvs,
            result_so_rsfc_la_mvs+result_so_rsfc_la_uvs,
            result_so_rsfc_ra_mvs+result_so_rsfc_ra_uvs,
            result_so_rsfc_lh_mvs+result_so_rsfc_lh_uvs,
            result_so_rsfc_rh_mvs+result_so_rsfc_rh_uvs,
            result_os_rsfc_la_mvs+result_os_rsfc_la_uvs,
            result_os_rsfc_ra_mvs+result_os_rsfc_ra_uvs,
            result_os_rsfc_lh_mvs+result_os_rsfc_lh_uvs,
            result_os_rsfc_rh_mvs+result_os_rsfc_rh_uvs
        )

        return results

    elif mode == 'rsFC-MMS':

        f_np = open('./bs_log/rsFC-MMS/NPP/results.pkl', 'rb')
        result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh, result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh, result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh, result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh = pickle.load(
            f_np)

        results_npp = (result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh, result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh,
                      result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh, result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh)

        f_p = open('./bs_log/rsFC-MMS/PP/results.pkl', 'rb')
        result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh, result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh, result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh, result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh = pickle.load(
            f_p)

        results_pp = (result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh, result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh,
                     result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh, result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh)

        f_pw = open('./bs_log/rsFC-MMS/PPw/results.pkl', 'rb')
        result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh, result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh, result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh, result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh = pickle.load(
            f_pw)

        results_ppw = (result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh, result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh,
                       result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh, result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh)

        f_psc = open('./bs_log/rsFC-MMS/CPP-SD/results.pkl', 'rb')
        result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh, result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh, result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh, result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh = pickle.load(
            f_psc)

        results_cpp_sd = (result_rsfc_la_mms_la, result_rsfc_la_mms_ra, result_rsfc_la_mms_lh, result_rsfc_la_mms_rh, result_rsfc_ra_mms_la, result_rsfc_ra_mms_ra, result_rsfc_ra_mms_lh, result_rsfc_ra_mms_rh,
                          result_rsfc_lh_mms_la, result_rsfc_lh_mms_ra, result_rsfc_lh_mms_lh, result_rsfc_lh_mms_rh, result_rsfc_rh_mms_la, result_rsfc_rh_mms_ra, result_rsfc_rh_mms_lh, result_rsfc_rh_mms_rh)

        return results_npp, results_pp, results_ppw, results_cpp_sd


def construct_parameter_configurations(mode):
    
    """
    Construct the whole parameter space for the specific combination

    Parameters
    ----------
    mode: string 'IMQ-MMS' or 'IMQ-rsFC' or 'rsFC-MMS'.
        The combination between two modalities.
    """

    if mode == 'IMQ-MMS':

        parameters_npp, parameters_pp, parameters_ppw, parameters_cpp_sd, combination, titles = [], [], [], [], [], []

        for dsf_x in DSF_X_mvs+DSF_X_uvs:
            for dsf_y in [
                    'pearson', 'euclidean', 'cosine', 'manhattan'
            ]:
                parameters_npp.append((dsf_x, dsf_y, None))

        for dsf_x in DSF_X_mvs+DSF_X_uvs:
            for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for P in [True]:
                    for PS in [(5, 5), (10, 10), (25, 25), (50, 50)]:
                        for PM in ['max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                            for SD in [None]:

                                mms = P, PS, PM, SD
                                parameters_pp.append((dsf_x, dsf_y, mms))
                        
        for dsf_x in DSF_X_mvs+DSF_X_uvs:
            for dsf_y in ['wrd', 'wmd']:
                for P in [True]:
                    for PS in [(25, 25), (50, 50)]:
                        for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                            for SD in [None]:

                                mms = P, PS, PM, SD
                                parameters_ppw.append((dsf_x, dsf_y, mms))
                        
        for dsf_x in ['pearson', 'euclidean', 'mahalanobis']:
            for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for P in [True]:
                    for PS in [(25, 25), (50, 50)]:
                        for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                            for SD in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:

                                mms = P, PS, PM, SD
                                parameters_cpp_sd.append((dsf_x, dsf_y, mms))
                        
        for dsf_x in ['cosine']:
            for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for P in [True]:
                    for PS in [(25, 25)]:
                        for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                            for SD in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:

                                mms = P, PS, PM, SD
                                parameters_cpp_sd.append((dsf_x, dsf_y, mms))
                        
        for dsf_x in ['cosine']:
            for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for P in [True]:
                    for PS in [(50, 50)]:
                        for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                            for SD in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:

                                mms = P, PS, PM, SD
                                parameters_cpp_sd.append((dsf_x, dsf_y, mms))
                        
        for dsf_x in ['manhattan']:
            for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for P in [True]:
                    for PS in [(25, 25), (50, 50)]:
                        for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                            for SD in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:

                                mms = P, PS, PM, SD
                                parameters_cpp_sd.append((dsf_x, dsf_y, mms))
                        
        for dsf_x in DSF_X_uvs:
            for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for P in [True]:
                    for PS in [(25, 25), (50, 50)]:
                        for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                            for SD in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:

                                mms = P, PS, PM, SD
                                parameters_cpp_sd.append((dsf_x, dsf_y, mms))

        for x in [(SS_mvs, SS_uvs), (SO_mvs, SO_uvs), (OS_mvs, OS_uvs)]:
            for y in [mms_la, mms_ra, mms_lh, mms_rh]:
                combination.append((x, y))

        for m in ['SS', 'SO', 'OS']:
            for mms in ['MMS-LA', 'MMS-RA', 'MMS-LH', 'MMS-RH']:

                titles.append(m + ' and '+mms)

        return parameters_npp, parameters_pp, parameters_ppw, parameters_cpp_sd, combination, titles

    elif mode == 'IMQ-rsFC':

        parameters, combination, titles = [], [], []

        for dsf_x in DSF_X_mvs+DSF_X_uvs:
            for dsf_y in [
                    'pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan'
            ]:
                parameters.append((dsf_x, dsf_y))

        for x in [(SS_mvs, SS_uvs), (SO_mvs, SO_uvs), (OS_mvs, OS_uvs)]:
            for y in [rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh]:
                combination.append((x, y))

        for m in ['SS', 'SO', 'OS']:
            for rsfc in ['rsFC-LA', 'rsFC-RA', 'rsFC-LH', 'rsFC-RH']:

                titles.append(m + ' and '+rsfc)

        return parameters, combination, titles

    elif mode == 'rsFC-MMS':

        parameters_npp, parameters_pp, parameters_ppw, parameters_cpp_sd, combination, titles = [], [], [], [], [], []

        for dsf_x in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
            for dsf_y in [
                    'pearson', 'euclidean', 'cosine', 'manhattan'
            ]:
                parameters_npp.append((dsf_x, dsf_y))
        
        for dsf_x in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
            for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for P in [True]:
                    for PS in [(5, 5), (10, 10), (25, 25), (50, 50)]:
                        for PM in ['max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                            for SD in [None]:

                                mms = P, PS, PM, SD
                                parameters_pp.append((dsf_x, dsf_y, mms))
                        
        for dsf_x in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
            for dsf_y in ['wrd', 'wmd']:
                for P in [True]:
                    for PS in [(25, 25), (50, 50)]:
                        for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                            for SD in [None]:
                        
                                mms = P, PS, PM, SD
                                parameters_ppw.append((dsf_x, dsf_y, mms))
                        
        for dsf_x in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
            for dsf_y in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                for P in [True]:
                    for PS in [(25, 25), (50, 50)]:
                        for PM in [None, 'max', 'mean', 'min', 'max-mean', 'max-min', 'mean-min', 'mmm']:
                            for SD in ['pearson', 'euclidean', 'mahalanobis', 'cosine', 'manhattan']:
                        
                                mms = P, PS, PM, SD
                                parameters_cpp_sd.append((dsf_x, dsf_y, mms))                       
        
        for x in [rsfc_la, rsfc_ra, rsfc_lh, rsfc_rh]:
            for y in [mms_la, mms_ra, mms_lh, mms_rh]:
                combination.append((x, y))
        
        for rsfc in ['rsFC-LA', 'rsFC-RA', 'rsFC-LH', 'rsFC-RH']:
            for mms in ['MMS-LA', 'MMS-RA', 'MMS-LH', 'MMS-RH']:

                titles.append(rsfc + ' and '+mms)

        return parameters_npp, parameters_pp, parameters_ppw, parameters_cpp_sd, combination, titles


def mean(result, n):

    """
    Compute the average IDM correspomndence from bootstrapping for each parameter configuration

    Parameters
    ----------
    result: list.
        The subject-wise bootstrapping results
    n: int. 
        The number of parameter configurations
    """

    return [np.mean(result[i][-1]) for i in range(n)]


def ks_test(data):

    """
    Perform two-sided Kolmogorov-Smirnov tests to test the normality of the distribution
    """
    
    print(kstest(np.array(data[0][-1])-np.array(data[1][-1]), cdf="norm"))
    print(kstest(np.array(data[0][-1])-np.array(data[2][-1]), cdf="norm"))
    print(kstest(np.array(data[0][-1])-np.array(data[3][-1]), cdf="norm"))
    print(kstest(np.array(data[1][-1])-np.array(data[2][-1]), cdf="norm"))
    print(kstest(np.array(data[1][-1])-np.array(data[3][-1]), cdf="norm"))
    print(kstest(np.array(data[2][-1])-np.array(data[3][-1]), cdf="norm"))
    

def wcx_test(data):

    """
    Perform two-sided Wilcoxon signed-rank tests to test for significant differences among brain regions
    """
    
    diffs, zvals, pvals = [], [], []
    
    stat_12 = np.mean(np.array(data[0][-1])-np.array(data[1][-1]))
    stat_13 = np.mean(np.array(data[0][-1])-np.array(data[2][-1]))
    stat_14 = np.mean(np.array(data[0][-1])-np.array(data[3][-1]))
    stat_23 = np.mean(np.array(data[1][-1])-np.array(data[2][-1]))
    stat_24 = np.mean(np.array(data[1][-1])-np.array(data[3][-1]))
    stat_34 = np.mean(np.array(data[2][-1])-np.array(data[3][-1]))
    
    l_12, h_12 = st.norm.interval(alpha=0.95,
                                  loc=stat_12,
                                  scale=st.sem(np.array(data[0][-1])-np.array(data[1][-1])))
    l_13, h_13 = st.norm.interval(alpha=0.95,
                                  loc=stat_13,
                                  scale=st.sem(np.array(data[0][-1])-np.array(data[2][-1])))
    l_14, h_14 = st.norm.interval(alpha=0.95,
                                  loc=stat_14,
                                  scale=st.sem(np.array(data[0][-1])-np.array(data[3][-1])))
    l_23, h_23 = st.norm.interval(alpha=0.95,
                                  loc=stat_23,
                                  scale=st.sem(np.array(data[1][-1])-np.array(data[2][-1])))
    l_24, h_24 = st.norm.interval(alpha=0.95,
                                  loc=stat_24,
                                  scale=st.sem(np.array(data[1][-1])-np.array(data[3][-1])))
    l_34, h_34 = st.norm.interval(alpha=0.95,
                                  loc=stat_34,
                                  scale=st.sem(np.array(data[2][-1])-np.array(data[3][-1])))
    
    r_12 = wilcoxon(data[0][-1], data[1][-1], alternative='two-sided', method='approx')
    r_13 = wilcoxon(data[0][-1], data[2][-1], alternative='two-sided', method='approx')
    r_14 = wilcoxon(data[0][-1], data[3][-1], alternative='two-sided', method='approx')
    r_23 = wilcoxon(data[1][-1], data[2][-1], alternative='two-sided', method='approx')
    r_24 = wilcoxon(data[1][-1], data[3][-1], alternative='two-sided', method='approx')
    r_34 = wilcoxon(data[2][-1], data[3][-1], alternative='two-sided', method='approx')
    
    diffs.append((stat_12, l_12, h_12)); diffs.append((stat_13, l_13, h_13)); diffs.append((stat_14, l_14, h_14)); diffs.append((stat_23, l_23, h_23)); diffs.append((stat_24, l_24, h_24)); diffs.append((stat_34, l_34, h_34))
    zvals.append(r_12.zstatistic); zvals.append(r_13.zstatistic); zvals.append(r_14.zstatistic); zvals.append(r_23.zstatistic); zvals.append(r_24.zstatistic); zvals.append(r_34.zstatistic)
    pvals.append(r_12.pvalue); pvals.append(r_13.pvalue); pvals.append(r_14.pvalue); pvals.append(r_23.pvalue); pvals.append(r_24.pvalue); pvals.append(r_34.pvalue)
        
    return diffs, zvals, pvals


def model_selection(results, combination, parameters, titles, mms_pc=None, show=False):
    
    """
    Select the optimal parameter configuration with the the highest average IDM correspondence from subject-wise bootstrapping
    
    Parameters
    ----------
    results: list.
        The subject-wise bootstrapping results grouped by the same pipeline.
    combination: list.
        The list of tuples, which consists of the combination of modalities.
    parameters: list.
        The list of tuples, which consists of the combination of parameters.
    titles: list.
        The list of tuples, which consists of the combination of names of modalities.
    mms_pc: string 'NPP' or 'PP' or 'PPw' or 'CPP-SD'. Default is None.
        The pipeline to construct inter-subject MMS dissimilarity matrix.
        'NPP' for without patching and pooling operations and computing surface distance,
        'PP' or 'PPw' (wmd and wrd) for with patching and pooling operations but without computing surface distance,
        'CPP-SD' for computing patching and pooling operations-based surface distance.
    bsm: bool True or False. Default is False.
        If bsm is True, then only return the highest average IDM correspondence
        If bsm is False, then return the info about the optimal parameter configuration
    show: bool True or False. Default is False.
        Print info or not.
    
    Returns
    -------
    winner: list. 
        The info about the optimal parameter configuration
    pvals: list. 
        The p-value about the optimal parameter configuration
    bs_mean: list. 
        The highest average IDM correspondence
    """

    winner, pvals, bs_mean = [], [], []
    num_p = len(parameters)

    for i, r in tqdm(enumerate(results)):

        idx, m = sorted([p for p in enumerate(mean(r, num_p))],
                        key=lambda x: x[-1],
                        reverse=True)[0]
        
        bs_mean.append((m, r[idx][-1]))
            
        x = combination[i][0]
                
        if 'SS' in titles[i] or 'SO' in titles[i] or 'OS' in titles[i]:
            
            if parameters[idx][0] in DSF_X_mvs:
                x = combination[i][0][0]
                
            if parameters[idx][0] in DSF_X_uvs:
                x = combination[i][0][1]
            
        if mms_pc == 'NPP':
            
            rho, pv = run(x,
                          combination[i][1],
                          parameters[idx][0],
                          parameters[idx][1],
                          iter=10000)

            if show:
                print(titles[i] + ': ' + 'rho ' + '%.4f' % float(rho) + ' p ' +
                      '%.4f' % float(pv) + ' achieved by ' + parameters[idx][0] +
                      '-' + parameters[idx][1] + ' ' + '%.4f' % float(m))

        else:
            rho, pv = run(x,
                          combination[i][1],
                          parameters[idx][0],
                          parameters[idx][1],
                          parameters[idx][2],
                          iter=10000)
                
            if show:
                print(titles[i] + ': ' + 'rho ' + '%.4f' % float(rho) + ' p ' +
                      '%.4f' % float(pv) + ' achieved by ' + parameters[idx][0] +
                      '-' + parameters[idx][1] + '-' +
                      str(parameters[idx][2][1][0]) + '-' +
                      str(parameters[idx][2][2]) + ' ' + '%.4f' % float(m))

        winner.append((r[idx], rho, m, parameters[idx]))
        pvals.append(pv)
            
    return winner, pvals, bs_mean


def eval_all(combination, w_npp, w_pp, w_ppw, w_cpp_sd, pv_npp, pv_pp, pv_ppw, pv_cpp_sd,
             titles):
    
    """
    Evaluate all the average IDM correspondence from bootstrapping for different pipelines
    """

    pvals_all = []

    for i in range(len(combination)):

        idx, m = sorted([(0, w_npp[i][2]), (1, w_pp[i][2]), (2, w_ppw[i][2]),
                         (3, w_cpp_sd[i][2])],
                        key=lambda x: x[-1],
                        reverse=True)[0]

        if idx == 0:
            rho = w_npp[i][1]
            pv = pv_npp[i]
            parameters = w_npp[i][-1]
            m = w_npp[i][2]
            l, h = st.norm.interval(alpha=0.95,
                                    loc=np.mean(w_npp[i][0][1]),
                                    scale=st.sem(w_npp[i][0][1]))
            print('NPP:    ' + titles[i] + ': ' + 'rho ' + '%.4f' % float(rho) +
                  ' p ' + '%.4f' % float(pv) + ' mean ' + '%.4f' % float(m) +
                  ' CI ' + '%.4f' % float(l) + '-' + '%.4f' % float(h) +
                  ' achieved by ' + parameters[0] + '-' + parameters[1])

        elif idx == 1:
            rho = w_pp[i][1]
            pv = pv_pp[i]
            parameters = w_pp[i][-1]
            m = w_pp[i][2]
            l, h = st.norm.interval(alpha=0.95,
                                    loc=np.mean(w_pp[i][0][1]),
                                    scale=st.sem(w_pp[i][0][1]))
            print('PP:     ' + titles[i] + ': ' + 'rho ' + '%.4f' % float(rho) +
                  ' p ' + '%.4f' % float(pv) + ' mean ' + '%.4f' % float(m) +
                  ' CI ' + '%.4f' % float(l) + '-' + '%.4f' % float(h) +
                  ' achieved by ' + parameters[0] + '-' + parameters[1] + '-' +
                  str(parameters[2][1][0]) + '-' + str(parameters[2][2]))

        elif idx == 2:
            rho = w_ppw[i][1]
            pv = pv_ppw[i]
            parameters = w_ppw[i][-1]
            m = w_ppw[i][2]
            l, h = st.norm.interval(alpha=0.95,
                                    loc=np.mean(w_ppw[i][0][1]),
                                    scale=st.sem(w_ppw[i][0][1]))
            print('PPw:    ' + titles[i] + ': ' + 'rho ' + '%.4f' % float(rho) +
                  ' p ' + '%.4f' % float(pv) + ' mean ' + '%.4f' % float(m) +
                  ' CI ' + '%.4f' % float(l) + '-' + '%.4f' % float(h) +
                  ' achieved by ' + parameters[0] + '-' + parameters[1] + '-' +
                  str(parameters[2][1][0]) + '-' + str(parameters[2][2]))

        elif idx == 3:
            rho = w_cpp_sd[i][1]
            pv = pv_cpp_sd[i]
            parameters = w_cpp_sd[i][-1]
            m = w_cpp_sd[i][2]
            l, h = st.norm.interval(alpha=0.95,
                                    loc=np.mean(w_cpp_sd[i][0][1]),
                                    scale=st.sem(w_cpp_sd[i][0][1]))
            print('CPP-SD: ' + titles[i] + ': ' + 'rho ' + '%.4f' % float(rho) +
                  ' p ' + '%.4f' % float(pv) + ' mean ' + '%.4f' % float(m) +
                  ' CI ' + '%.4f' % float(l) + '-' + '%.4f' % float(h) +
                  ' achieved by ' + parameters[0] + '-' + parameters[1] + '-' +
                  str(parameters[2][1][0]) + '-' + str(parameters[2][2]) + '-' + str(parameters[2][3]))

        pvals_all.append(pv)

    return pvals_all


def data_prep(results_mms, results_rsfc, target, combination_mms, combination_rsfc, parameters_mms, parameters_rsfc, titles):

    """
    Prepare the dataset for the following dyadic regression analysis
    """
    
    imq_idms, mms_idms, rsfc_idms = [], [], []
    num_p_mms = len(parameters_mms)
    num_p_rsfc = len(parameters_rsfc)

    print('IMQ-MMS')

    for i, r in enumerate(results_mms):
        
        r = [(j[0], [-1]) if j[0].split('-')[0] != target else j for j in r]

        idx, m = sorted([p for p in enumerate(mean(r, num_p_mms))],
                        key=lambda x: x[-1],
                        reverse=True)[0]
            
        if parameters_mms[idx][0] in DSF_X_mvs:
            x = combination_mms[i][0][0]
                
        if parameters_mms[idx][0] in DSF_X_uvs:
            x = combination_mms[i][0][1]
            
        if parameters_mms[idx][-1] == None:            
            rho, imq_idm, mms_idm = run(x,
                                        combination_mms[i][1],
                                        parameters_mms[idx][0],
                                        parameters_mms[idx][1],
                                        idm=True)

            print(titles[i] + ': ' + '%.4f' % float(rho) + ' achieved by ' + parameters_mms[idx][0] +
                  '-' + parameters_mms[idx][1] + ' ' + '%.4f' % float(m))

        else:
            rho, imq_idm, mms_idm = run(x,
                                        combination_mms[i][1],
                                        parameters_mms[idx][0],
                                        parameters_mms[idx][1],
                                        parameters_mms[idx][2],
                                        idm=True)

            print(titles[i] + ': ' + '%.4f' % float(rho) + ' achieved by ' + parameters_mms[idx][0] +
                  '-' + parameters_mms[idx][1] + '-' +
                  str(parameters_mms[idx][2][1][0]) + '-' +
                  str(parameters_mms[idx][2][2]) + ' ' + '%.4f' % float(m))

        imq_idms.append(imq_idm)
        mms_idms.append(mms_idm)
    
    print('------------------------------------------------')
    print('IMQ-rsFC')
        
    for i, r in enumerate(results_rsfc):
        
        r = [(j[0], [-1]) if j[0].split('-')[0] != target else j for j in r]

        idx, m = sorted([p for p in enumerate(mean(r, num_p_rsfc))],
                        key=lambda x: x[-1],
                        reverse=True)[0]
            
        x = combination_rsfc[i][0]
            
        if parameters_rsfc[idx][0] in DSF_X_mvs:
            x = combination_rsfc[i][0][0]
                
        if parameters_rsfc[idx][0] in DSF_X_uvs:
            x = combination_rsfc[i][0][1]
            
        rho, _, rsfc_idm = run(x,
                               combination_rsfc[i][1],
                               parameters_rsfc[idx][0],
                               parameters_rsfc[idx][1],
                               idm=True)

        print(titles[i] + ': ' + '%.4f' % float(rho) + ' achieved by ' + parameters_rsfc[idx][0] +
              '-' + parameters_rsfc[idx][1] + ' ' + '%.4f' % float(m))

        rsfc_idms.append(rsfc_idm)
    
    subject_left = [i+1 for i in range(24) for j in range(24) if i < j]+[j+1 for i in range(24) for j in range(24) if i < j]
    subject_right = [j+1 for i in range(24) for j in range(24) if i < j]+[i+1 for i in range(24) for j in range(24) if i < j]
    
    for i, t in enumerate(titles):
    
        df = DataFrame({
            t[:2]: list(zscore(imq_idms[i], nan_policy='omit')*-1)+list(zscore(imq_idms[i], nan_policy='omit')*-1),
            t[-2:]+'_MMS': list(zscore(mms_idms[i])*-1)+list(zscore(mms_idms[i])*-1),
            t[-2:]+'_rsFC': list(zscore(rsfc_idms[i])*-1)+list(zscore(rsfc_idms[i])*-1),
            'Subject_left': subject_left,
            'Subject_right': subject_right,
        })
        df.to_csv('./data/dra/'+t+'/'+target+'.csv')


def regression_model(dv, iv):
    
    """
    Construct all possible regerssion models (i.e., iterate all the distance metrics)

    Parameters
    ----------
    dv: string 'SS' or 'SO' or 'OS'.
        'SS' for self-self mentalisation, 'SO' for self-other mentalisation, 'OS' for other-self mentalisation.
    iv: string 'RH' or 'LA' or 'RA'.
        'RH' for right hippocampus, 'LA' for left amygdala, 'RA' for right amygdala.
    """

    for dsf in DSF_X_mvs + DSF_X_uvs:

        data = pd.read_csv('./data/dra/' + dv + '-' +
                           iv + '/' + dsf + '.csv')
        print(
            '----------------------------------------------------------------------------------------------------'
        )
        print(dsf)
        print(spearmanr(data[iv + "_MMS"], data[iv + "_rsFC"]))
        try:
            model = Lmer(dv + " ~ " + iv + "_MMS * " + iv +
                         "_rsFC + (1|Subject_left) + (1|Subject_right)",
                         data=data)
            print(model.fit())
        except:
            pass


def specific_model(dsf, dv, iv):
    
    """
    Construct the optimal regression model, which has the highest log-likelihood
    
    Parameters
    ----------
    dsf: string 'abs' or 'mean' or 'min' or 'abs*mean' or 'max' or 'pearson' or 'euclidean' or 'mahalanobis' or 'cosine' or 'manhattan'.
        The method to calculate mentalising similarity.
    dv: string 'SS' or 'SO' or 'OS'.
        'SS' for self-self mentalisation, 'SO' for self-other mentalisation, 'OS' for other-self mentalisation.
    iv: string 'RH' or 'LA' or 'RA'.
        'RH' for right hippocampus, 'LA' for left amygdala, 'RA' for right amygdala.
    """

    data = pd.read_csv('./data/dra/' + dv + '-' + iv +
                       '/' + dsf + '.csv')

    model = Lmer(dv + " ~ " + iv + "_MMS * " + iv +
                 "_rsFC + (1|Subject_left) + (1|Subject_right)",
                 data=data)

    min_rsFC = np.quantile(data[iv + '_rsFC'], 0.25)
    max_rsFC = np.quantile(data[iv + '_rsFC'], 0.75)

    categorical_rsFC = []

    for i in data[iv + '_rsFC'][:276]:

        if np.abs(i - min_rsFC) < np.abs(i - max_rsFC):
            categorical_rsFC.append('Q1')
        else:
            categorical_rsFC.append('Q3')

    print(model.fit())     
    
    return model, list(data[iv + '_MMS'][:276]), list(data[iv + '_rsFC'][:276]), list(model.fits[:276]), categorical_rsFC