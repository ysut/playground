import numpy as np
import pandas as pd


class GenePrioritization:
    def __init__(self, args: dict):
        self._root_path = args['resources']
        self.path_to_gado = f'{self._root_path}/GADO_score.tsv'

    def __str__(self):
        return f'{self.gene} {self.score} {self.rank}'



# Input a GADO result
