import pandas as pd
import glob

# root_path = '/work-Github/Developing/WES/db/gnomAD/'

class Gnomad:
    def __init__(self, root_path):
        self._root_path = root_path
        self.path_to_gnomad = glob.glob(f'{self._root_path}/gnomAD/*.txt') # path_to_gnomad is list

    def annotation(self, dataframe):
        gnomad = pd.read_table(self.path_to_gnomad[0], header=0, dtype=str)
        gnomad = gnomad[['gene', 'pLI', 'pRec', 'syn_z', 'mis_z', 'oe_lof_upper']]
        gnomad = gnomad.rename(columns={'gene': 'gnomAD_gene', 'oe_lof_upper': 'LOEUF'})
        dataframe = pd.merge(dataframe, gnomad, left_on='Gene.refGene', right_on='gnomAD_gene', how='left')

        return dataframe
