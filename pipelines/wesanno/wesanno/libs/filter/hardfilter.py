import pandas as pd
from .gtfilter import ModelDataFrame


class HardFilter:
    """Hard filtration
    This class is to filter out possibly non-pathogenic variants.

    1. Exclude unmatched inheritance mode
       'MOI' column shows well-known inheritance mode from HGMD and G2P.
        
    2. Unreported genes filtration
       'DM' column shows the number of variants registerd as DM class in HGMD.
       The variants on unreported genes (DM is 0) are filterd out.
       
    3. Exclude the variants called with possibly low quality
       'GQ' and 'DP' columns show the quality of variants.
       Then, variants with low quality are filtered out.

    """
    def __init__(self, dfs: ModelDataFrame):
        self.dfs = dfs


    def __exclude_known_moi(self, df: pd.DataFrame, mode: str) -> pd.DataFrame:
        if ((mode == 'XL') |(mode == 'YL')):
            df.loc[((df['MOI'] == 'XL') 
                    | (df['MOI'] == 'YL') 
                    | (df['MOI'] == '.')),
                    'MOI_FILTER'] = 'PASS'      
        else:
            df.loc[((df['MOI'] == mode) 
                    | (df['MOI'] == '.')),
                    'MOI_FILTER'] = 'PASS'
        
        return df


    def hard_filtering(self) -> ModelDataFrame:
        self.dfs.AD = self.__exclude_known_moi(self.dfs.AD, mode='AD')
        self.dfs.Hm = self.__exclude_known_moi(self.dfs.Hm, mode='AR')
        self.dfs.CH = self.__exclude_known_moi(self.dfs.CH, mode='AR')
        self.dfs.XL = self.__exclude_known_moi(self.dfs.XL, mode='XL')

        for df in [self.dfs.AD, self.dfs.Hm, self.dfs.CH, self.dfs.XL]:
            df.loc[
                (df['MOI_FILTER'] == 'PASS') 
                & (df['non_zero_DM'] == 'PASS')
                & (df['GQ_FILTER'] == 'PASS')
                # & (df['DP_FILTER'] == 'PASS')
                & (df['FlaggedSNP_FILTER'] == 'PASS'),
                'HARD_FILTER'] = 'PASS'
                    
        return self.dfs
