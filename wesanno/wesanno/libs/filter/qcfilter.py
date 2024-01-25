import pandas as pd

class QcFilter:
    def __init__(self, df):
        self.df = df

    def __exclude_low_gq(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[
            (df['GQ'] < 20),
            'GQ_FILTER'
            ] = 'PASS'

        return df
    
    def __exclude_low_dp(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[
            (df['DP'] < 2),
            'DP_FILTER'
            ] = 'PASS'

        return df
    
    def __exclude_low_ab(self, df: pd.DataFrame) -> pd.DataFrame:
        df['AB'] = df['AD'] / df['DP']
        
        df.loc[
            (df['AB'] < 0.2),
            'AB_FILTER'
            ] = 'PASS'

        return df
    
    def __exclude_low_ad(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[
            (df['AD'] < 5),
            'AD_FILTER'
            ] = 'PASS'

        return df
    
    def exclude_low_qc(self) -> pd.DataFrame:
        self.df = self.__exclude_low_gq(self.df)
        self.df = self.__exclude_low_dp(self.df)

        return self.df