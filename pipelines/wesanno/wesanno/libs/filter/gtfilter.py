import sys
from dataclasses import dataclass

import pandas as pd

@dataclass
class ModelDataFrame():
    AD: pd.DataFrame
    Hm: pd.DataFrame
    CH: pd.DataFrame
    XL: pd.DataFrame

class GtFilter:
    def __init__(self, df: pd.DataFrame, mode_samples_info):
        self.df = df
        self.mode = mode_samples_info.mode


    def __separate_chrom(self, df: pd.DataFrame) -> tuple:
        dfA = df[(df['CHROM'] != 'X') & (df['CHROM'] != 'Y')].copy()
        dfS = df[(df['CHROM'] == 'X') | (df['CHROM'] == 'Y')].copy()

        return dfA, dfS


    def __hetero_count(self, df: pd.DataFrame) -> pd.DataFrame:
        sr = df.groupby('Gene.refGene')['Gene.refGene'].transform('count').rename('Hetero_count')
        df = pd.concat([df, sr], axis=1)
        
        return df


    def _gt_ad_filter(self, df: pd.DataFrame):
        if self.mode == 'proband':
            df = self.__hetero_count(df)
            df.loc[
                (df['GT_Pro'] == '0/1') 
                & (df['Hetero_count'] == 1),
                'GT_FILTER'
                ] = 'PASS'

            return df

        elif self.mode == 'duo':
            df = self.__hetero_count(df)
            df.loc[
                (~df['GT_Pro'].str.contains('0/0|\./\.')) 
                & (df['GT_Par'].str.contains('0/0|\./\.')) 
                & (df['Hetero_count'] == 1),
                'GT_FILTER'
                ] = 'PASS'

            return df
        
        elif self.mode == 'trio':
            df.loc[
                (~df['GT_Pro'].str.contains('0/0|\./\.')) 
                & (df['GT_Fa'].str.contains('0/0|\./\.')) 
                & (df['GT_Mo'].str.contains('0/0|\./\.')),
                'GT_FILTER'
                ] = 'PASS'

            return df

        elif self.mode == 'quad_affected':
            df.loc[
                (~df['GT_Pro'].str.contains('0/0|\./\.'))
                & (df['GT_Fa'].str.contains('0/0|\./\.'))
                & (df['GT_Mo'].str.contains('0/0|\./\.'))
                & ~(df['GT_Sib'].str.contains('0/0|\./\.')),
                'GT_FILTER'
                ] = 'PASS'

            return df

        elif self.mode == 'quad_unaffected':
            df.loc[
                (~df['GT_Pro'].str.contains('0/0|\./\.'))
                & (df['GT_Fa'].str.contains('0/0|\./\.'))
                & (df['GT_Mo'].str.contains('0/0|\./\.'))
                & (df['GT_Sib'].str.contains('0/0|\./\.')),
                'GT_FILTER'
                ] = 'PASS'

            return df
        
        else:
            print('Error: Mode is not correct.')
            sys.exit(1)


    def __gt_hm_filter(self, df):
        if ((self.mode == 'proband') or (self.mode == 'duo')):
            df.loc[
                (df['GT_Pro'] == '1/1'),
                'GT_FILTER'
                ] = 'PASS'
            
            return df

        elif self.mode == 'trio':
            df. loc[
                (df['GT_Pro'] == '1/1') 
                & (df['GT_Fa'] != '1/1') 
                & (df['GT_Mo'] != '1/1'),
                'GT_FILTER'
                ] = 'PASS'

            return df
        
        elif self.mode == 'quad_affected':
            df.loc[
                (df['GT_Pro'] == '1/1') 
                & (df['GT_Fa'] != '1/1') 
                & (df['GT_Mo'] != '1/1')
                & (df['GT_Sib'] == '1/1'),
                'GT_FILTER'
                ] = 'PASS'
        
        elif self.mode == 'quad_unaffected':
            df.loc[
                (df['GT_Pro'] == '1/1') 
                & (df['GT_Fa'] != '1/1') 
                & (df['GT_Mo'] != '1/1')
                & (df['GT_Sib'] != '1/1'),
                'GT_FILTER'
                ] = 'PASS'

            return df
        
        else:
            print('Error: Mode is not correct.')
            sys.exit(1)


    def __gt_ch_filter(self, df):
        df = self.__hetero_count(df)
        df.loc[(df['GT_Pro'] == '0/1') & (df['Hetero_count'] >= 2),
                'GT_FILTER'] = 'PASS'

        return df


    def __gt_xl_filter(self, df) -> pd.DataFrame:
        df.loc[~df['GT_Pro'].str.contains('0/0|\./\.'), 
                'GT_FILTER'] = 'PASS'
        
        return df


    def genotypeing_filter(self) -> ModelDataFrame:
        dfA, dfS = self.__separate_chrom(self.df)

        dfs = ModelDataFrame(
            AD=self._gt_ad_filter(dfA),
            Hm=self.__gt_hm_filter(dfA),
            CH=self.__gt_ch_filter(dfA),
            XL=self.__gt_xl_filter(dfS)
        )

        return dfs