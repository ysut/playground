from dataclasses import dataclass
import glob

import numpy as np
import pandas as pd
import sqlite3

# Set up logger
from logging import getLogger, StreamHandler, Formatter, INFO

logger = getLogger(__name__)
format = Formatter(fmt='%(asctime)s [%(levelname)s] - %(message)s',
                   datefmt='%Y-%m-%d %H:%M:%S')
handler = StreamHandler()
handler.setLevel(INFO)
handler.setFormatter(format)
logger.addHandler(handler)
logger.setLevel(INFO)
logger.propagate = False


@dataclass
class DFs:
    df: pd.DataFrame
    df_old: pd.DataFrame
    df_mailed: pd.DataFrame


def find_excel_files() -> tuple:
    logger.info("Find the excel files")
    old_excel: str = glob.glob('rawdata/2004*.xlsx')[0]
    new_excel: str = glob.glob('rawdata/自動報告書*')[0]
    mailed_excel: str = glob.glob('rawdata/郵送*')[0]

    return old_excel, new_excel, mailed_excel


def load_excel_file(old_excel: str, new_excel: str, mailed_excel: str) -> DFs:
    ## Columns to use
    usecols_new = ['Date', 'YCU ID', 'DNA ID', '病名', 'State', 'Gene', 
                   'WES batch', '家族関係', 'father DNA ID', 'mother DNA ID', 
                   'Proband_DNA_ID']

    usecols_old = ['Date', 'YCU ID', 'DNA ID', '病名', 'WES batch', 
                   '家族関係', 'father DNA ID', 'mother DNA ID', 
                   'Proband_DNA_ID']

    usecols_post = ['検体拝受日', 'YCUID', 'DNAID', '診断名', 'State', 
                    '解析結果(遺伝子名)', 'ProbandID']

    # Load excel files as pd.DataFrame
    logger.info(f"Load {new_excel} file")
    df = pd.read_excel(new_excel, dtype=str, sheet_name='database', 
                    index_col=None, skiprows=3, usecols=usecols_new)
    df.dropna(how='all')
    df = df[usecols_new]

    logger.info(f"Load {old_excel} file")
    df_old = pd.read_excel(old_excel, index_col=None, usecols=usecols_old)
    df_old.dropna(how='all')
    df_old = df_old[usecols_old]

    logger.info(f"Load {mailed_excel} file")
    df_mailed = pd.read_excel(mailed_excel, index_col=None, usecols=usecols_post)
    df_mailed.dropna(how='all')
    df_mailed = df_mailed[usecols_post]

    rename_dict = {
        'YCU ID': 'YCU_ID', 'YCUID': 'YCU_ID', 
        'DNA ID': 'DNA_ID', 'DNAID': 'DNA_ID'
        }
    df.rename(columns=rename_dict, inplace=True)
    df_old.rename(columns=rename_dict, inplace=True)
    df_mailed.rename(columns=rename_dict, inplace=True)

    return DFs(df, df_old, df_mailed)


def adjust_state_info(df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Adjust 'State' information")
    replace_to_identified = [
        '89717692', '22006383', '166210771', '18602433', '18598085', '53279465', 
        '6495539', '89653823', '89712015', '37840968', '47990943', '24024725', 
        '39929311', '29,610,505', '50607726', '50607674', '75959244', 
        '220154743', '220156581', '56385389', '8242896', '166909430', 
        '153363092', '166237229', '62070967', '105834449', '130422388', 
        '1736004', '52082841', '101953473', '101953473', 'Confirm', 
        'SNV_identified', 'CNV_identified', 'SV_identified', 
        'SV_determined', 'identified', 'Undetermined→identified', 
        'Repeat_expansion', 'repeat expansion identified', 'mosicism susp',
        'SNV_identified;CNV_identified', 'Tandem_repeat_expansion_identified'
        ]

    replace_to_undetermined = [
        'On going', 'On_going', 'undetermined', '検体取り下げ', 
        '欠番', 'Unknown', 'Inconclusive', 
        'サンガーの結果、本家系で見られたCOL4A2バリアントはなかった'
        ]
    
    df.fillna({'State': 'Undetermined'}, inplace=True)
    df.replace({'State': replace_to_identified}, 'Identified', inplace=True)
    df.replace({'State': replace_to_undetermined}, 'Undetermined', inplace=True)

    return df


def create_db(db_path: str, dfs: DFs) -> None:
    # Connect to sqlite3 database
    logger.info("Connect to sqlite3 database")
    conn = sqlite3.connect(db_path)

    # Create sqlite database (If exists, replace). Three tables will be created.
    logger.info("Create a sqlite database in the 'db' directory")
    dfs.df.to_sql('new_samples', conn, if_exists='replace', index=False)
    dfs.df_old.to_sql('old_samples', conn, if_exists='replace', index=False)
    dfs.df_mailed.to_sql('mailed_samples', conn, if_exists='replace', index=False)

    conn.close()

    return None


def main():
    # Loading excel files as pd.DataFrame
    old_excel, new_excel, mailed_excel = find_excel_files()
    dfs = load_excel_file(old_excel, new_excel, mailed_excel)

    # Adjust 'State' information
    dfs.df = adjust_state_info(dfs.df)
    dfs.df_mailed = adjust_state_info(dfs.df_mailed)

    # Fill NaN with np.nan
    logger.info("Fill NaN with np.nan")
    for df in [dfs.df, dfs.df_old, dfs.df_mailed]:
        dfs.df.replace('', np.nan, inplace=True)

    if len(dfs.df['State'].unique()) > 2:
        errmsg = (f"State column has more than 2 unique values. "
                  f"In the {new_excel}.")
        logger.error(errmsg)
        logger.error(dfs.df['State'].unique())
        exit(1)

    if len(dfs.df_mailed['State'].unique()) > 2:
        errmsg = (f"State column has more than 2 unique values. "
                  f"In the {mailed_excel}.")
        logger.error(errmsg)
        logger.error(dfs.df_mailed['State'].unique())
        exit(1)

    # Create a db
    ycu_db = 'db/ycudb.db'
    create_db(ycu_db, dfs)



if __name__ == '__main__':
    main()
    logger.info("Done")
