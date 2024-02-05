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
            ((df['ExAC_ALL'] < 0.001) | (df['ExAC_ALL'].isnull()))
            & ((df['ExAC_EAS'] < 0.001) | (df['ExAC_EAS'].isnull()))
            & ((df['ExAC_FIN'] < 0.001) | (df['ExAC_FIN'].isnull()))
            & ((df['ExAC_NFE'] < 0.001) | (df['ExAC_NFE'].isnull()))
            & ((df['ExAC_OTH'] < 0.001) | (df['ExAC_OTH'].isnull()))
            & ((df['ExAC_SAS'] < 0.001) | (df['ExAC_SAS'].isnull()))
            & ((df['ExAC_AFR'] < 0.001) | (df['ExAC_AFR'].isnull()))
            & ((df['ExAC_AMR'] < 0.001) | (df['ExAC_AMR'].isnull()))
            & ((df['HGVD_AF'] < 0.001) | (df['HGVD_AF'].isnull()))
            & ((df['ABraOM_AF'] < 0.001) | (df['ABraOM_AF'].isnull()))
            & ((df['esp6500siv2_all'] < 0.001) | (df['esp6500siv2_all'].isnull()))
            & ((df['ToMMo3.5KJPN_AF'] < 0.001) | (df['ToMMo3.5KJPN_AF'].isnull())),
            'MAF_0.1%_FILTER'
            ] = 'PASS'

        return self.df


    def __extract_1_percent(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[
            ((df['ExAC_ALL'] < 0.01) | (df['ExAC_ALL'].isnull()))
            & ((df['ExAC_EAS'] < 0.01) | (df['ExAC_EAS'].isnull()))
            & ((df['ExAC_FIN'] < 0.01) | (df['ExAC_FIN'].isnull()))
            & ((df['ExAC_NFE'] < 0.01) | (df['ExAC_NFE'].isnull()))
            & ((df['ExAC_OTH'] < 0.01) | (df['ExAC_OTH'].isnull()))
            & ((df['ExAC_SAS'] < 0.01) | (df['ExAC_SAS'].isnull()))
            & ((df['ExAC_AFR'] < 0.01) | (df['ExAC_AFR'].isnull()))
            & ((df['ExAC_AMR'] < 0.01) | (df['ExAC_AMR'].isnull()))
            & ((df['HGVD_AF'] < 0.01) | (df['HGVD_AF'].isnull()))
            & ((df['ABraOM_AF'] < 0.01) | (df['ABraOM_AF'].isnull()))
            & ((df['esp6500siv2_all'] < 0.01) | (df['esp6500siv2_all'].isnull()))
            & ((df['ToMMo3.5KJPN_AF'] < 0.01) | (df['ToMMo3.5KJPN_AF'].isnull())),
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
        self.df = self.__extract_zero_inhouse(self.df)
        self.df = self.__extract_1_percent_inhouse(self.df)
        self.df = self.__extract_01_percent(self.df)
        self.df = self.__extract_1_percent(self.df)
        self.df = self.__exclude_nonflagged_snps(self.df)

        return self.df
