import pandas as pd

from logging import getLogger
logger = getLogger(__name__)

class QcFilter:
    def __init__(self, df: pd.DataFrame, configs: dict):
        self.gq: int = configs['Quality_cutoffs']['GQ']
        self.dp: int = configs['Quality_cutoffs']['DP']
        self.df: pd.DataFrame = df
        logger.info(f"FILTER cutoffs: GQ >= {self.gq}  DP >= {self.dp}")

    def __exclude_low_gq(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[df['GQ'] >= self.gq, 'GQ_FILTER'] = 'PASS'

        return df    

    def __exclude_low_dp(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[df['DP'] >= self.dp, 'DP_FILTER'] = 'PASS'

        return df
        
    def exclude_low_quality(self) -> pd.DataFrame:
        self.df = self.__exclude_low_gq(self.df)
        # self.df = self.__exclude_low_dp(self.df)

        return self.df