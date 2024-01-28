import dataclasses
import pandas as pd
import numpy as np
from liftover import get_lifter

def liftover_to_hg38(dfs: dataclasses.dataclass, anno_sheets: list
                     ) -> dataclasses.dataclass:
    converter = get_lifter('hg19', 'hg38')
    for sheet in anno_sheets:
        dfs.sheets[sheet]['POS_hg38'] = converter.convert_coordinate(
            dfs.sheets[sheet]['CHROM'], dfs.sheets[sheet]['POS'])
        
def split_alt(dfs: dataclasses.dataclass, anno_sheets: list
              ) -> dataclasses.dataclass:
    pass