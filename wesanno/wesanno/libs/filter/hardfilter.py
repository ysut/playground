import pandas as pd
from .gtfilter import ModelDataFrame


class HardFilter:
    def __init__(self, dfs: ModelDataFrame):
        self.dfs = dfs
        

    def __exclude_known_moi(self, df: pd.DataFrame, mode: str) -> pd.DataFrame:

        if ((mode == 'XL') |(mode == 'YL')):
            df = df[((df['MOI'] == 'XL') 
                     | (df['MOI'] == 'YL') 
                     | (df['MOI'] == '.'))]        
        else:
            df = df[((df['MOI'] == mode) 
                     | (df['MOI'] == '.'))]
        
        return df

    def __exclude_low_gqdp(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df[df['GQ_FILTER'] == 'PASS']
        df = df[df['DP_FILTER'] == 'PASS']

        return df

    def __exclude_dm_zero(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df[df['DM'] != 0]

        return df

    def hard_filtering(self) -> ModelDataFrame:
        self.dfs.AD = self.__exclude_known_moi(self.dfs.AD, mode='AD')
        self.dfs.Hm = self.__exclude_known_moi(self.dfs.Hm, mode='AR')
        self.dfs.CH = self.__exclude_known_moi(self.dfs.CH, mode='AR')
        self.dfs.XL = self.__exclude_known_moi(self.dfs.XL, mode='XL')

        self.dfs.AD = self.__exclude_low_gqdp(self.dfs.AD)
        self.dfs.Hm = self.__exclude_low_gqdp(self.dfs.Hm)
        self.dfs.CH = self.__exclude_low_gqdp(self.dfs.CH)
        self.dfs.XL = self.__exclude_low_gqdp(self.dfs.XL)

        self.dfs.AD = self.__exclude_dm_zero(self.dfs.AD)
        self.dfs.Hm = self.__exclude_dm_zero(self.dfs.Hm)
        self.dfs.CH = self.__exclude_dm_zero(self.dfs.CH)
        self.dfs.XL = self.__exclude_dm_zero(self.dfs.XL)
                
        return self.dfs
