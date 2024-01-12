import pandas as pd
import glob

# root_path = '/work-Github/Developing/WES/db/DECIPHER/'

class Decipher:
    def __init__(self, root_path):
        self._root_path = root_path
        self.path_to_ddg2p = glob.glob(f'{self._root_path}/DECIPHER/DDG2P_*.csv')
        self.path_to_eyeg2p = glob.glob(f'{self._root_path}/DECIPHER/EyeG2P_*.csv')
        self.path_to_sking2p = glob.glob(f'{self._root_path}/DECIPHER/SkinG2P_*.csv')
        self.path_to_cardiacg2p = glob.glob(f'{self._root_path}/DECIPHER/CardiacG2P_*.csv')

    def annotation(self, dataframe):
        ddg2p = pd.read_csv(self.path_to_ddg2p[0], header=0, dtype=str)
        eyeg2p = pd.read_csv(self.path_to_eyeg2p[0], header=0, dtype=str)
        sking2p = pd.read_csv(self.path_to_sking2p[0], header=0, dtype=str)
        cardiacg2p = pd.read_csv(self.path_to_cardiacg2p[0], header=0, dtype=str)

        dcprDropList = ['gene mim', 'disease name', 'disease mim', 'panel',
                        'phenotypes', 'organ specificity list', 'pmids', 
                        'prev symbols', 'hgnc id', 'gene disease pair entry date', 
                        'mutation consequence flag', 'confidence value flag', 
                        'comments', 'variant consequence', 'disease ontology']
        dbs = ['DDG2P', 'EyeG2P', 'SkinG2P', 'CardiacG2P']
        dfs = [ddg2p, eyeg2p, sking2p, cardiacg2p]

        for (df, db) in zip(dfs, dbs):
            df = df.drop(columns=dcprDropList).add_prefix(f'{db}_')
            dataframe = pd.merge(dataframe, df, how='left',
                                 left_on='Gene.refGene', right_on=f'{db}_gene symbol')
        dataframe.drop(columns=['DDG2P_gene symbol', 'EyeG2P_gene symbol',
                                'SkinG2P_gene symbol', 'CardiacG2P_gene symbol'], inplace=True)

        return dataframe
