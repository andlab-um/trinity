# -*- coding: utf-8 -*-

__author__ = 'Zhaoning Li'

from util import *
import scipy.stats as st
from statsmodels.stats.multitest import multipletests
from scipy.stats import kstest, wilcoxon

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
    
    diffs, pvals = [], []
    
    stat_12 = np.mean(np.array(data[0][-1])-np.array(data[1][-1]))
    stat_13 = np.mean(np.array(data[0][-1])-np.array(data[2][-1]))
    stat_14 = np.mean(np.array(data[0][-1])-np.array(data[3][-1]))
    stat_23 = np.mean(np.array(data[1][-1])-np.array(data[2][-1]))
    stat_24 = np.mean(np.array(data[1][-1])-np.array(data[3][-1]))
    stat_34 = np.mean(np.array(data[2][-1])-np.array(data[3][-1]))
    
    _, pval_12 = wilcoxon(data[0][-1], data[1][-1], alternative='two-sided')
    _, pval_13 = wilcoxon(data[0][-1], data[2][-1], alternative='two-sided')
    _, pval_14 = wilcoxon(data[0][-1], data[3][-1], alternative='two-sided')
    _, pval_23 = wilcoxon(data[1][-1], data[2][-1], alternative='two-sided')
    _, pval_24 = wilcoxon(data[1][-1], data[3][-1], alternative='two-sided')
    _, pval_34 = wilcoxon(data[2][-1], data[3][-1], alternative='two-sided')
    
    diffs.append(stat_12); diffs.append(stat_13); diffs.append(stat_14); diffs.append(stat_23); diffs.append(stat_24); diffs.append(stat_34) 
    pvals.append(pval_12); pvals.append(pval_13); pvals.append(pval_14); pvals.append(pval_23); pvals.append(pval_24); pvals.append(pval_34)
        
    return diffs, pvals









