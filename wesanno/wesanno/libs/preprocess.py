import os
import re
import sys
import numpy as np
import pandas as pd
from typing import NamedTuple
from pandarallel import pandarallel

# from liftover import get_lifter

from logging import getLogger
logger = getLogger(__name__)

os.environ['JOBLIB_TEMP_FOLDER'] = '/tmp' 
pandarallel.initialize(progress_bar=True, verbose=2)

class PreProcessExomeSummary:
    def __init__(self, df: pd.DataFrame, args: dict, mode_samples_info):
        self.df = df
        self.args = args
        self.mode = mode_samples_info.mode
        self.mode_samples_info = mode_samples_info

    def __extract_hgvd_maf(
            self, df: pd.DataFrame, col='HGVD') -> pd.DataFrame:
        df['HGVD_AF'] = df[col].str.extract('(\.|[0-1]\.\d{,6})')
        df['HGVD_AF'] = df['HGVD_AF'].replace('.', np.nan)
        df['HGVD_AF'] = df['HGVD_AF'].astype(float)
        return df
    
    def __extract_tommo_maf(
            self, df: pd.DataFrame, 
            col='snp20171005_tommo3.5k_passed') -> pd.DataFrame:
        df['ToMMo3.5KJPN_AF'] = df[col].str.extract('(\.|[0-1]\.\d{0,})')
        df['ToMMo3.5KJPN_AF'] = df['ToMMo3.5KJPN_AF'].replace('.', np.nan)
        df['ToMMo3.5KJPN_AF'] = df['ToMMo3.5KJPN_AF'].astype(float)
        return df

    def __extract_abraom_maf(
        self, df: pd.DataFrame, col='ABRaOM') -> pd.DataFrame:
        if col in df.columns:
            df['ABraOM_AF'] = df[col].str.extract('(?<=Frequencies=)(\d\.\d+)')
            df['ABraOM_AF'] = df['ABraOM_AF'].astype(float)
        else:
            df['ABraOM_AF'] = np.nan

        return df
    
    def __extract_inhouse_maf_auto(self, row) -> str:
        if ((row['CHROM'] != 'X') & (row['CHROM'] != 'Y')):
            return re.search(r'(?<=Mutation=.:575:)(\d+)', row['INFO']).group()
        else:
            return '.'
        
    def __extract_inhouse_maf_m(self, row) -> str:
        if ((row['CHROM'] == 'X') | (row['CHROM'] == 'Y')):
            mf = re.search(r'(?<=Mutation=.:281,294:)(.+)', row['INFO']).group()
            return mf.split(',')[0]
        else:
            return '.'

    def __extract_inhouse_maf_f(self, row) -> str:
        if ((row['CHROM'] == 'X') | (row['CHROM'] == 'Y')):
            mf = re.search(r'(?<=Mutation=.:281,294:)(.+)', row['INFO']).group()
            return mf.split(',')[1]
        else:
            return '.'
    
    def __extract_genotyoeing_info(
            self, df: pd.DataFrame) -> pd.DataFrame:
        if self.mode == 'proband':
            proband_id = self.mode_samples_info.proband_id
            df['GT_Pro'] = df[proband_id].str.extract('(\./\.|\d/\d)')
        elif self.mode == 'duo':
            proband_id = self.mode_samples_info.proband_id
            parent_id = self.mode_samples_info.parent_id
            df['GT_Pro'] = df[proband_id].str.extract('(\./\.|\d/\d)')
            df['GT_Par'] = df[parent_id].str.extract('(\./\.|\d/\d)')
        elif self.mode == 'trio':
            proband_id = self.mode_samples_info.proband_id
            father_id = self.mode_samples_info.father_id
            mother_id = self.mode_samples_info.mother_id
            df['GT_Pro'] = df[proband_id].str.extract('(\./\.|\d/\d)')
            df['GT_Fa'] = df[father_id].str.extract('(\./\.|\d/\d)')
            df['GT_Mo'] = df[mother_id].str.extract('(\./\.|\d/\d)')
        elif ((self.mode == 'quad_affected') 
              | (self.mode == 'quad_unaffected')):
            proband_id = self.mode_samples_info.proband_id
            father_id = self.mode_samples_info.father_id
            mother_id = self.mode_samples_info.mother_id
            sibling_id = self.mode_samples_info.sibling_id
            df['GT_Pro'] = df[proband_id].str.extract('(\./\.|\d/\d)')
            df['GT_Fa'] = df[father_id].str.extract('(\./\.|\d/\d)')
            df['GT_Mo'] = df[mother_id].str.extract('(\./\.|\d/\d)') 
            df['GT_Sib'] = df[sibling_id].str.extract('(\./\.|\d/\d)') 
        else:
            print('Erorr: Mode is not correct.')
            sys.exit(1)
        
        return df


    def __replace_dot_to_nan(self, df: pd.DataFrame) -> pd.DataFrame:
        df['ExAC_ALL'].replace('.', np.nan, inplace=True)
        df['ExAC_EAS'].replace('.', np.nan, inplace=True)
        df['ExAC_NFE'].replace('.', np.nan, inplace=True)
        df['ExAC_SAS'].replace('.', np.nan, inplace=True)
        df['ExAC_AMR'].replace('.', np.nan, inplace=True)
        df['ExAC_AFR'].replace('.', np.nan, inplace=True)
        df['ExAC_FIN'].replace('.', np.nan, inplace=True)
        df['ExAC_OTH'].replace('.', np.nan, inplace=True)
        df['esp6500siv2_all'].replace('.', np.nan, inplace=True)
        df['InHouse575_AF'].replace('.', np.nan, inplace=True)
        df['InHouseMale_AF'].replace('.', np.nan, inplace=True)
        df['InHouseFemale_AF'].replace('.', np.nan, inplace=True)
        df = df.astype(
            {'ExAC_ALL': float, 'ExAC_EAS': float, 'ExAC_NFE': float, 
             'ExAC_SAS': float, 'ExAC_AMR': float, 'ExAC_AFR': float, 
             'ExAC_FIN': float, 'ExAC_OTH': float, 'esp6500siv2_all': float,
             'InHouse575_AF': float, 'InHouseMale_AF': float, 
             'InHouseFemale_AF': float}
             )
        
        return df

    def _split_qc_info(self, df: pd.DataFrame) -> pd.DataFrame:
        proband_id = self.mode_samples_info.proband_id
        for i, new_col in enumerate(['GT', 'AD', 'DP', 'GQ', 'PL']):
            df[new_col] = df[proband_id].str.split(':').str[i]
        
        df.replace({'GQ': {'.': 0}}, inplace=True)
        df.fillna({'GQ': 0}, inplace=True)
        df = df.astype({'GQ': 'int32'})

        df.replace({'DP': {'.': 0}}, inplace=True)
        df.fillna({'DP': 0}, inplace=True)
        df = df.astype({'DP': 'int32'})
    
        return df
    
    def __process_dataframe_optimized(self, df: pd.DataFrame) -> pd.DataFrame:
        # Save the original column names
        column_names: list = list(df.columns) + ['SplitALT']

        # Initialize a list to store the processing results
        rows_list = []

        # Iterate over the rows of the dataframe (Maybe low performance ......)
        for row in df.itertuples(index=False):
            # Split the value of the REF column by comma
            alts = row.ALT.split(',')
            # Generate a new row for each split value
            for alt in alts:
                new_row = row._asdict()
                new_row['SplitALT'] = alt
                rows_list.append(new_row)

        # Convert the list to a dataframe
        result_df = pd.DataFrame(rows_list)
        result_df.columns = column_names 

        return result_df

    def __drop_unused_cols(self, df: pd.DataFrame) -> pd.DataFrame:
        droplist = [
            'exac03_pLI_Z:pLI:pRec', 'dbscSNV_RF_SCORE', 'dpsi_max_tissue', 
            'dpsi_zscore', 'Otherinfo1', 'OtherInfo2', 'OtherInfo3',
            'HGVD', 'snp20171005_tommo3.5k_passed', 'dbscSNV_ADA_SCORE',
            'gnomAD_v2_1_1_June2020_z_pLI:LOUEF:pRec', 'GT'
            ]
        
        for dropcolumn in droplist:
            if dropcolumn in df.columns:
                pass
            else:
                logger.info(f"{dropcolumn} is not in the dataframe.")
                droplist.remove(dropcolumn)

        return df.drop(droplist, axis=1)
    

    ### Public methods ###
    def liftover_to_hg38(self, row) -> int:
        converter = get_lifter('hg19', 'hg38')
        result = converter.query(row['CHROM'], int(row['POS']))
        if result:
            return int(result[0][1])
        else:
            return None


    def all_pre_processing(self):
        # if ((self.args['assembly'] == 'hg19') 
        #     | (self.args['assembly'] == 'GRCh37')):
        #     logger.info('Liftover to hg38 is started ...')
        #     logger.info('It takes about 5 minutes.')
        #     self.df['POS_hg38'] = self.df.parallel_apply(
        #         self.liftover_to_hg38, axis=1)
        #     logger.info('Liftover to hg38 is finished.')
        # else:
        #     logger.info(f"No liftover to hg38 "
        #                 f"(Assembly is {self.args['assembly']}).")
        #     pass

        # Extract InHouse MAF
        logger.info('Extract InHouse MAF')
        self.df['InHouse575_AF'] = self.df.parallel_apply(
            self.__extract_inhouse_maf_auto, axis=1)
        self.df['InHouseMale_AF'] = self.df.parallel_apply(
            self.__extract_inhouse_maf_m, axis=1)
        self.df['InHouseFemale_AF'] = self.df.parallel_apply(
            self.__extract_inhouse_maf_f, axis=1)

        # Extract MAF from each database cols
        logger.info('Extract MAF from each database cols')
        self.df = self.__extract_hgvd_maf(self.df)
        self.df = self.__extract_tommo_maf(self.df)
        self.df = self.__extract_abraom_maf(self.df)

        # Extract genotypeing info
        logger.info('Extract genotypeing info')
        self.df = self.__extract_genotyoeing_info(self.df)

        # Split QC info
        logger.info('Split QC info')
        self.df = self._split_qc_info(self.df)

        # Replace '.' to np.nan
        logger.info('Replace "." to np.nan')
        self.df = self.__replace_dot_to_nan(self.df)

        # Split ALT column
        logger.info('Split ALT column')
        self.df = self.__process_dataframe_optimized(self.df)
        self.df = self.df.drop_duplicates(
            subset=['CHROM', 'POS', 'ALT', 'SplitALT']
            )

        # Drop unused columns
        logger.info('Drop unused columns')
        self.df = self.__drop_unused_cols(self.df)

        return self.df