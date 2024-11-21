import pandas as pd
from typing import NamedTuple

class MafFilter:
    def __init__(self, df, mode_samples_info, configs):
        self.df = df
        self.mode = mode_samples_info.mode
        self.config = configs


    def __extract_1_percent_inhouse(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[
            ((df['InHouse575_AF'] < 6 ) | (df['InHouse575_AF'].isnull()))
            & ((df['InHouseMale_AF'] < 3) | (df['InHouseMale_AF'].isnull()))
            & ((df['InHouseFemale_AF'] < 3) | (df['InHouseFemale_AF'].isnull())),
            'InHouse_1%_FILTER'
            ] = 'PASS'

        return self.df


    def __extract_zero_inhouse(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[
            ((df['InHouse575_AF'] == 0 ) | (df['InHouse575_AF'].isnull()))
            & ((df['InHouseMale_AF'] == 0) | (df['InHouseMale_AF'].isnull()))
            & ((df['InHouseFemale_AF'] == 0) | (df['InHouseFemale_AF'].isnull())),
            'InHouse_absent_FILTER'
            ] = 'PASS'

        return self.df


    def __extract_01_percent(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[
            ((df['PopFreqMax'] < 0.001) | (df['PopFreqMax'].isnull()))
            & ((df['gnomADv4_exome_AF_popmax'] < 0.001) | (df['gnomADv4_exome_AF_popmax'].isnull()))
            & ((df['gnomADv4_genome_AF_popmax'] < 0.001) | (df['gnomADv4_genome_AF_popmax'].isnull()))
            & ((df['HGVD_AF'] < 0.001) | (df['HGVD_AF'].isnull()))
            & ((df['ABraOM_AF'] < 0.001) | (df['ABraOM_AF'].isnull()))
            & ((df['ToMMo54K'] < 0.001) | (df['ToMMo54K'].isnull()))
            & ((df['ToMMo54K_PAR2'] < 0.001) | (df['ToMMo54K_PAR2'].isnull()))
            & ((df['ToMMo54K_PAR3'] < 0.001) | (df['ToMMo54K_PAR3'].isnull()))
            & ((df['ToMMo8.3K_MT_MALE'] < 0.001) | (df['ToMMo8.3K_MT_MALE'].isnull()))
            & ((df['ToMMo8.3K_MT_FEMALE'] < 0.001) | (df['ToMMo8.3K_MT_FEMALE'].isnull())),      
            'MAF_0.1%_FILTER'
            ] = 'PASS'

        return self.df


    def __extract_1_percent(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[
            ((df['PopFreqMax'] < 0.01) | (df['PopFreqMax'].isnull()))
            & ((df['gnomADv4_exome_AF_popmax'] < 0.01) | (df['gnomADv4_exome_AF_popmax'].isnull()))
            & ((df['gnomADv4_genome_AF_popmax'] < 0.01) | (df['gnomADv4_genome_AF_popmax'].isnull()))
            & ((df['HGVD_AF'] < 0.01) | (df['HGVD_AF'].isnull()))
            & ((df['ABraOM_AF'] < 0.01) | (df['ABraOM_AF'].isnull()))
            & ((df['ToMMo54K'] < 0.01) | (df['ToMMo54K'].isnull()))
            & ((df['ToMMo54K_PAR2'] < 0.01) | (df['ToMMo54K_PAR2'].isnull()))
            & ((df['ToMMo54K_PAR3'] < 0.01) | (df['ToMMo54K_PAR3'].isnull()))
            & ((df['ToMMo8.3K_MT_MALE'] < 0.01) | (df['ToMMo8.3K_MT_MALE'].isnull()))
            & ((df['ToMMo8.3K_MT_FEMALE'] < 0.01) | (df['ToMMo8.3K_MT_FEMALE'].isnull())),  
            'MAF_1%_FILTER'
            ] = 'PASS'

        return self.df


    def __extract_custom_maf_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        """
        # cutoffs = self.conifgs
        return df


    def __exclude_nonflagged_snps(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[df['snp138NonFlagged'] == '.', 'FlaggedSNP_FILTER'] = 'PASS'
        
        return df


    def all_filtering(self):
        # self.df = self.__extract_zero_inhouse(self.df)
        # self.df = self.__extract_1_percent_inhouse(self.df)
        self.df = self.__extract_01_percent(self.df)
        self.df = self.__extract_1_percent(self.df)
        self.df = self.__exclude_nonflagged_snps(self.df)

        return self.df
