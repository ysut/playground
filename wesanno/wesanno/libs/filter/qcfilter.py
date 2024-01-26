import pandas as pd

from logging import getLogger
logger = getLogger(__name__)

class QcFilter:
    def __init__(self, configs: dict):
        self.gq: int = configs['Quality_cutoffs']['GQ']
        self.dp: int = configs['Quality_cutoffs']['DP']
        logger.info(f"FILTER cutoffs: GQ >= {self.gq}  DP >= {self.dp}")


    def __exclude_low_gq(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[df['GQ'] >= self.gq, 'GQ_FILTER'] = 'PASS'

        return df
    
    def __exclude_low_dp(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[df['DP'] >= self.dp, 'DP_FILTER'] = 'PASS'

        return df
    
    def __exclude_low_ab(self, df: pd.DataFrame) -> pd.DataFrame: 
        # Allele balance filter
        # AB>=0.2 for heterozygotes
        pass
    
    
    def exclude_low_quality(self, df:pd.DataFrame) -> pd.DataFrame:
        df = self.__exclude_low_gq(df)
        df = self.__exclude_low_dp(df)

        return df