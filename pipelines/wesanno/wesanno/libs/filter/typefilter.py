import pandas as pd

class TypeFilter:
    def __init__(self, df):
        self.df = df


    def __exclude_hlamuc(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[
            ~df['Gene.refGene'].str.contains('HLA|MUC'),
            'HLAMUC_FILTER'
            ] = 'PASS'

        return df


    def __exclude_exonic_synonymous(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[
            ~((df['ExonicFunc.refGene'] == 'synonymous SNV') 
            & (df['Func.refGene'] == 'exonic')),
            'ExonicSyno_FILTER'
            ] = 'PASS'

        return df


    def exclude_hlamuc_and_exonicsyno(self) -> pd.DataFrame:
        self.df = self.__exclude_hlamuc(self.df)
        self.df = self.__exclude_exonic_synonymous(self.df)

        return self.df