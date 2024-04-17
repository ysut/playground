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
pandarallel.initialize(progress_bar=False, verbose=1)

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

    def __extract_tommo_maf(self, df: pd.DataFrame):
        df['ToMMo54K'] = df['INFO'].str.extract('(?<=ToMMo54KJPN_lifted37_AF=)(\d\.?\d*)')
        df['ToMMo54K_PAR2'] = df['INFO'].str.extract('(?<=ToMMo54KJPN_lifted37_AF_PAR2=)(\d\.?\d*)')
        df['ToMMo54K_PAR3'] = df['INFO'].str.extract('(?<=ToMMo54KJPN_lifted37_AF_PAR3=)(\d\.?\d*)')
        df['ToMMo8.3K_MT_MALE'] = df['INFO'].str.extract('(?<=ToMMo8.3KJPN_AF_MALE_MT=)(\d\.?\d*)')
        df['ToMMo8.3K_MT_FEMALE'] = df['INFO'].str.extract('(?<=ToMMo8.3KJPN_AF_FEMALE_MT=)(\d\.?\d*)')
        return df
    
    # def __extract_tommo_maf(
    #         self, df: pd.DataFrame, 
    #         col='snp20171005_tommo3.5k_passed') -> pd.DataFrame:
    #     df['ToMMo3.5KJPN_AF'] = df[col].str.extract('(\.|[0-1]\.\d{0,})')
    #     df['ToMMo3.5KJPN_AF'] = df['ToMMo3.5KJPN_AF'].replace('.', np.nan)
    #     df['ToMMo3.5KJPN_AF'] = df['ToMMo3.5KJPN_AF'].astype(float)
    #     return df

    def __extract_spliceai(self, df: pd.DataFrame):
        df['SpliceAI_raw'] = df['INFO'].str.extract('(?<=SpliceAI=)([^;]+)')
        df['SpliceAI_tmp'] = df['SpliceAI_raw'].str.split(',').str[0]
        df = pd.concat([df, df['SpliceAI_tmp'].str.split('|', expand=True)], axis=1)
        df.drop(columns=['SpliceAI_tmp', 0], inplace=True)
        df.rename(columns={1: 'SpliceAI_Tx', 
                        2: 'DS_AG', 3: 'DS_AL', 4: 'DS_DG', 5: 'DS_DL', 
                        6: 'DP_AG', 7: 'DP_AL', 8: 'DP_DG', 9: 'DP_DL'}, 
                        inplace=True)
        df['SpliceAI_Max'] = df[['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']].max(axis=1)
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
            mf = re.search(r'(?<=Mutation=.:281,294:)([^;]+)', row['INFO']).group()
            return mf.split(',')[0]
        else:
            return '.'

    def __extract_inhouse_maf_f(self, row) -> str:
        if ((row['CHROM'] == 'X') | (row['CHROM'] == 'Y')):
            mf = re.search(r'(?<=Mutation=.:281,294:)([^;]+)', row['INFO']).group()
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

    def __rename_maf_cols(self, df: pd.DataFrame) -> pd.DataFrame:
        
        gnomad_cols = [
            'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 
            'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin','AF_asj', 'AF_oth',
            'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 
            'non_cancer_AF_popmax', 'controls_AF_popmax']
        
        for col in df.columns:
            if col in gnomad_cols:
                new_col = 'gnomADv4_exome_' + col
                df = df.rename(columns={col: new_col})
            elif col.lstrip('2_') in gnomad_cols:
                new_col = 'gnomADv4_genome_' + col.lstrip('2_')
                df = df.rename(columns={col: new_col})
            else:
                pass
        
        return df

    def __replace_dot_to_nan(self, df: pd.DataFrame) -> pd.DataFrame:
        maf_cols = [
            'gnomADv4_exome_AF', 'gnomADv4_exome_AF_popmax', 
            'gnomADv4_exome_AF_male', 'gnomADv4_exome_AF_female', 
            'gnomADv4_exome_AF_raw', 'gnomADv4_exome_AF_afr', 
            'gnomADv4_exome_AF_sas', 'gnomADv4_exome_AF_amr', 
            'gnomADv4_exome_AF_eas', 'gnomADv4_exome_AF_nfe', 
            'gnomADv4_exome_AF_fin', 'gnomADv4_exome_AF_asj', 
            'gnomADv4_exome_AF_oth', 
            'gnomADv4_exome_non_topmed_AF_popmax',
            'gnomADv4_exome_non_neuro_AF_popmax', 
            'gnomADv4_exome_non_cancer_AF_popmax',
            'gnomADv4_exome_controls_AF_popmax', 
            'gnomADv4_genome_AF', 'gnomADv4_genome_AF_popmax', 
            'gnomADv4_genome_AF_male', 'gnomADv4_genome_AF_female', 
            'gnomADv4_genome_AF_raw', 'gnomADv4_genome_AF_afr', 
            'gnomADv4_genome_AF_amr', 'gnomADv4_genome_AF_eas', 
            'gnomADv4_genome_AF_nfe', 'gnomADv4_genome_AF_fin', 
            'gnomADv4_genome_AF_asj', 'gnomADv4_genome_AF_oth', 
            'gnomADv4_genome_non_topmed_AF_popmax',
            'gnomADv4_genome_non_neuro_AF_popmax', 
            'gnomADv4_genome_controls_AF_popmax', 
            'PopFreqMax', 'ExAC_ALL', 'ESP6500siv2_ALL', 'CG46', 
            'GME_AF', 'ABraOM_AF', 
            'InHouse575_AF', 'InHouseMale_AF', 'InHouseFemale_AF',
            'HGVD_AF', 'ToMMo54K', 'ToMMo54K_PAR2', 'ToMMo54K_PAR3',
            'ToMMo8.3K_MT_MALE', 'ToMMo8.3K_MT_FEMALE']

        for col in maf_cols:
            df[col] = df[col].replace('.', np.nan)
            df[col] = df[col].fillna(np.nan)
            df[col] = df[col].astype(float)
        
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
    
    def __split_alt_col(self, df: pd.DataFrame) -> pd.DataFrame:
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
        unused_cols = [
            'exac03_pLI_Z:pLI:pRec', 'dbscSNV_RF_SCORE', 'dpsi_max_tissue', 
            'dpsi_zscore', 'Otherinfo1', 'Otherinfo2', 'Otherinfo3', 
            'OtherInfo2', 'OtherInfo3', 'HGVD', 
            'snp20171005_tommo3.5k_passed', 'dbscSNV_ADA_SCORE',
            'gnomAD_v2_1_1_June2020_z_pLI:LOUEF:pRec', 'GT',
            'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE',
            'ExAC_OTH', 'ExAC_SAS', 'ESP6500siv2_AA', 'ESP6500siv2_EA', 
            '1000G_AFR', '1000G_AMR', '1000G_EAS', '1000G_EUR', '1000G_SAS',
            'GME_NWA', 'GME_NEA', 'GME_AP', 'GME_Israel', 'GME_SD', 'GME_TP',
            'GME_CA', 'gnomADv4_exomeAF_raw', 'gnomADv4_genome_AF_raw', 
            'Polyphen2_HDIV_score', 'Polyphen2_HDIV_rankscore', 
            'Polyphen2_HDIV_pred', 
            ]
        droplist = []

        for col in unused_cols:
            if col not in df.columns:
                pass
            else:
                droplist.append(col)

        return df.drop(droplist, axis=1)
    

    ### Public methods ###
    def liftover_to_hg38(self, row) -> int:
        converter = get_lifter('hg19', 'hg38')
        result = converter.query(row['CHROM'], int(row['POS']))
        if result:
            return int(result[0][1])
        else:
            return None


    ### Main method ###
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

        # Drop ALL dot rows columns
        logger.info('Drop ALL dot rows columns')
        self.df = self.df.loc[:,~(self.df == '.').all(axis=0)]

        # Extract InHouse MAF
        logger.info('Extract InHouse MAF')
        self.df.loc[:,'InHouse575_AF'] = self.df.parallel_apply(
            self.__extract_inhouse_maf_auto, axis=1)
        self.df.loc[:,'InHouseMale_AF'] = self.df.parallel_apply(
            self.__extract_inhouse_maf_m, axis=1)
        self.df.loc[:,'InHouseFemale_AF'] = self.df.parallel_apply(
            self.__extract_inhouse_maf_f, axis=1)

        # Extract MAF from each database cols
        logger.info('Extract MAF from each database cols')
        self.df = self.__extract_hgvd_maf(self.df).copy()
        self.df = self.__extract_tommo_maf(self.df).copy()
        self.df = self.__extract_abraom_maf(self.df)

        # Extract genotypeing info
        logger.info('Extract genotypeing info')
        self.df = self.__extract_genotyoeing_info(self.df)

        # Split QC info
        logger.info('Split QC info')
        self.df = self._split_qc_info(self.df)

        # Rename gnomAD MAF cols
        logger.info('Rename gnomAD MAF cols')
        self.df = self.__rename_maf_cols(self.df)

        # Replace '.' to np.nan
        logger.info('Replace "." to np.nan')
        self.df = self.__replace_dot_to_nan(self.df)

        # Split ALT column
        logger.info('Split ALT column')
        self.df = self.__split_alt_col(self.df)
        self.df = self.df.drop_duplicates(
            subset=['CHROM', 'POS', 'ALT', 'SplitALT']
            )
        
        # Extract SpliceAI
        logger.info('Extract SpliceAI')
        self.df = self.__extract_spliceai(self.df)

        # Drop unused columns
        logger.info('Drop unused columns')
        self.df = self.__drop_unused_cols(self.df)

        return self.df