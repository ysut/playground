import os

import pandas as pd
from pandarallel import pandarallel
import numpy as np
from liftover import ChainFile
from logging import getLogger

from hyperlink import ExcelSheetsDF

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
    # Chain file path in the Docker container is hard-coded
    chain_file = '/app/resources/hg19ToHg38.over.chain.gz'
    converter = ChainFile(chain_file, 'hg19', 'hg38')
    try:   
        result = converter[row['CHROM']][row['POS']][0][1]
        return result
    except IndexError:
        logger.warning(f"Failed to liftover ({row['CHROM']}:{row['POS']})")
        return np.nan


def liftover_to_hg38(
        dfs:ExcelSheetsDF, args: dict, anno_sheets: list) ->ExcelSheetsDF:
    for sheet in anno_sheets:
        df = dfs.sheets[sheet]
        if is_liftover_needed(df, args):
            logger.info(f"Inserting liftover columns ('POS_hg38') to {sheet}")
            df['POS_hg38'] = df.parallel_apply(liftover_process, axis=1)
        else:
            logger.info(f"Skip liftover columns insertion to {sheet}")
    
    return dfs


def is_split_alt_needed(args: dict) -> bool:
    if args['alt_col'] == 'ALT':
        logger.info('Skip split ALT column step. Using "ALT" column as is.')
        return False
    else:
        if args['splitALTFlag']:
            return True
        else:
            logger.info('Skip split ALT column step. Using specified ALT column as is.')
            return False
        

def split_alt_process(dfs:ExcelSheetsDF, anno_sheets: list) ->ExcelSheetsDF:
    for sheet in anno_sheets:
        df = dfs.sheets[sheet]
        columns: list = list(df.columns) + ['SplitALT']
        rows_list = []
        for row in df.itertuples(index=False):
            alts = row.ALT.split(',')
            for alt in alts:
                new_row = row._asdict()
                new_row['ALT'] = alt
                rows_list.append(new_row)
        new_df = pd.DataFrame(rows_list)
        new_df.columns = columns
        dfs.sheets[sheet] = new_df

    return dfs


def split_alt_col(
        dfs:ExcelSheetsDF, args: dict, anno_sheets: list) ->ExcelSheetsDF:
    if is_split_alt_needed(args):
        logger.info('Split ALT column processing ......')
        dfs = split_alt_process(dfs, anno_sheets)
    else:
        pass
    
    return dfs
