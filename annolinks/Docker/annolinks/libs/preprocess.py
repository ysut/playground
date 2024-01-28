import dataclasses
import os

import pandas as pd
from pandarallel import pandarallel
import numpy as np
from liftover import get_lifter
from logging import getLogger

pandarallel.initialize(
    progress_bar=True, verbose=2, use_memory_fs=False,
    nb_workers=int(os.cpu_count() -1))
os.environ['JOBLIB_TEMP_FOLDER'] = '/tmp' 
logger = getLogger(__name__)

def is_liftover_needed(df: pd.DataFrame, args: dict) -> bool:
    if 'DECIPHER' in args['skip_sites']:
        return False
    else:
        if ((args['assembly'] == 'hg38') | (args['assembly'] == 'GRCh38')):
            return False
        else:
            if args['pos38'] in df.columns:
                return False
            else:
                return True
            
            
def liftover_process(row) -> int:
    converter = get_lifter('hg19', 'hg38')
    return converter[row['CHROM']][row['POS']][0][1]


def liftover_to_hg38(dfs: dataclasses.dataclass, args: dict,
                      anno_sheets: list) -> dataclasses.dataclass:
    for sheet in anno_sheets:
        df = dfs.sheets[sheet]
        if is_liftover_needed(df, args):
            logger.info(f"Inserting liftover columns ('POS_hg38') to {sheet}")
            df['POS_hg38'] = df.parallel_apply(liftover_process, axis=1)
        else:
            logger.info(f"Skip liftover columns insertion to {sheet}")
    
    return dfs


def split_alt(dfs: dataclasses.dataclass, anno_sheets: list
              ) -> dataclasses.dataclass:
    pass