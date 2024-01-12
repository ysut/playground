import pandas as pd
import glob

# root_path = '/work-Github/Developing/WES/db/gnomAD/'

class Hgmd:
    def __init__(self, root_path):
        self._root_path = root_path
        self.path_to_hgmd = glob.glob(f'{self._root_path}/HGMD/*.txt') # path_to_gnomad is list

    def annotation(self, dataframe):
        hgmd = pd.read_table(self.path_to_hgmd[0], header=0, dtype=str)
        hgmd = hgmd[['GeneSymbol', 'expected_inheritance', 'DM']]
        hgmd = hgmd.rename(columns={'expected_inheritance': 'exp.MOI', 'DM': 'DM'})
        dataframe = pd.merge(dataframe, hgmd, left_on='gnomAD_gene', right_on='GeneSymbol', how='left')

        return dataframe
