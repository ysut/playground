import pandas as pd
from .gtfilter import ModelDataFrame


from logging import getLogger
logger = getLogger(__name__)
print(__name__)

def counter(dfs: ModelDataFrame, output_excel: str) -> ModelDataFrame:
    # Pick up columns to count variants
    # cols_for_counter = [
    #     'InHouse_absent_FILTER', 'InHouse_1%_FILTER', 
    #     'MAF_0.1%_FILTER', 'MAF_1%_FILTER', 
    #     'HLAMUC_FILTER', 'ExonicSyno_FILTER', 'GT_FILTER'
    #     ]
    # dfs.AD = dfs.AD[cols_for_counter]
    # dfs.Hm = dfs.Hm[cols_for_counter]
    # dfs.CH = dfs.CH[cols_for_counter]
    # dfs.XL = dfs.XL[cols_for_counter]
    
    # Count all variants
    num_all_vars: int = len(dfs.AD) + len(dfs.XL)
    all_vars: dict = {
        'AD': num_all_vars,
        'Homo': num_all_vars,
        'CH': num_all_vars,
        'XL': num_all_vars
        }

    # CHROM FILTER
    post_chrom: dict = {
        'AD': len(dfs.AD), 
        'Homo': len(dfs.Hm), 
        'CH': len(dfs.CH), 
        'XL': len(dfs.XL)
        }

    # GT FILTER
    dfs.AD = dfs.AD[dfs.AD['GT_FILTER'] == 'PASS']
    dfs.Hm = dfs.Hm[dfs.Hm['GT_FILTER'] == 'PASS']
    dfs.CH = dfs.CH[dfs.CH['GT_FILTER'] == 'PASS']
    dfs.XL = dfs.XL[dfs.XL['GT_FILTER'] == 'PASS']
    post_gt: dict = {
        'AD': len(dfs.AD), 
        'Homo': len(dfs.Hm), 
        'CH': len(dfs.CH), 
        'XL': len(dfs.XL)
        }

    # MAF FILTER
    dfs.AD = dfs.AD[dfs.AD['MAF_0.1%_FILTER'] == 'PASS']
    dfs.Hm = dfs.Hm[dfs.Hm['MAF_1%_FILTER'] == 'PASS']
    dfs.CH = dfs.CH[dfs.CH['MAF_1%_FILTER'] == 'PASS']
    dfs.XL = dfs.XL[dfs.XL['MAF_1%_FILTER'] == 'PASS']
    post_maf: dict = {
        'AD': len(dfs.AD), 
        'Homo': len(dfs.Hm), 
        'CH': len(dfs.CH), 
        'XL': len(dfs.XL)
        }

    # InHouse FILTER
    # dfs.AD = dfs.AD[dfs.AD['InHouse_absent_FILTER'] == 'PASS']
    # dfs.Hm = dfs.Hm[dfs.Hm['InHouse_1%_FILTER'] == 'PASS']
    # dfs.CH = dfs.CH[dfs.CH['InHouse_1%_FILTER'] == 'PASS']
    # dfs.XL = dfs.XL[dfs.XL['InHouse_1%_FILTER'] == 'PASS']
    post_inhouse: dict = {
        'AD': len(dfs.AD), 
        'Homo': len(dfs.Hm), 
        'CH': len(dfs.CH), 
        'XL': len(dfs.XL)
        }

    # HLA/MUC FILTER
    dfs.AD = dfs.AD[dfs.AD['HLAMUC_FILTER'] == 'PASS']
    dfs.Hm = dfs.Hm[dfs.Hm['HLAMUC_FILTER'] == 'PASS']
    dfs.CH = dfs.CH[dfs.CH['HLAMUC_FILTER'] == 'PASS']
    dfs.XL = dfs.XL[dfs.XL['HLAMUC_FILTER'] == 'PASS']
    post_hlamuc: dict = {
        'AD': len(dfs.AD), 
        'Homo': len(dfs.Hm), 
        'CH': len(dfs.CH), 
        'XL': len(dfs.XL)
        }

    # ExonicSyno FILTER
    dfs.AD = dfs.AD[dfs.AD['ExonicSyno_FILTER'] == 'PASS']
    dfs.Hm = dfs.Hm[dfs.Hm['ExonicSyno_FILTER'] == 'PASS']
    dfs.CH = dfs.CH[dfs.CH['ExonicSyno_FILTER'] == 'PASS']
    dfs.XL = dfs.XL[dfs.XL['ExonicSyno_FILTER'] == 'PASS']
    post_exsyno: dict = {
        'AD': len(dfs.AD), 
        'Homo': len(dfs.Hm), 
        'CH': len(dfs.CH), 
        'XL': len(dfs.XL)
        }
    
    # Generate dataframe of counter result
    index: list = [
        'All variants',
        'Chrom. Filter', 
        'GT Fitler', 
        'MAF Filter', 
        'In-House Filter',
        'Exclude HLA/MUC',
        'Exclude Ex.Syno.'
        ]
    
    counter_summary = pd.DataFrame(
        [all_vars, post_chrom, post_gt, post_maf, 
        post_inhouse, post_hlamuc, post_exsyno], index=index)

    counter_summary.to_excel(output_excel, sheet_name='FilteringSummary')

    return dfs

