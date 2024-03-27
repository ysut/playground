import dataclasses
import sys
import numpy as np
import pandas as pd
import toml
import pathlib2 as pathlib

with open("config.toml") as f:
    config = toml.load(f)

@dataclasses.dataclass
class Thresholds:
    GQ: int = config['QC_Thresholds']['GQ']
    DP: int = config['QC_Thresholds']['DP']
    AB_het = config['QC_Thresholds']['AB_het']
    AB_hom = config['QC_Thresholds']['AB_hom']
    SIFT: float = config['In_silico_predictions']['SIFT']
    PP2: float = config['In_silico_predictions']['PolyPhen-2']
    CADD: float = config['In_silico_predictions']['CADD']
    MutationTaster: list = dataclasses.field(default_factory=list)
    GGM_AF: float = config['AF_Thresholds']['GGM']

# Define the thresholds of 'MutationTaster' as a list of strings
Thresholds.MutationTaster = config['In_silico_predictions']['exclude_MutationTaster']


def _split_qc_col(df: pd.DataFrame) -> pd.DataFrame:
    for rel in ['pro', 'pat', 'mat']:
        df = pd.concat(
            [df, df[f'GQ:DP:AD({rel})'].str.split(':', expand=True)], axis=1)
        for i in range(3):
            df[i] = df[i].replace('.', np.nan)
            df[i] = df[i].replace('-', np.nan)
        df = df.astype({0: float, 1: float, 2: float})
        df.rename(columns={0: f'GQ({rel})', 1: f'DP({rel})', 2: f'AD({rel})'}, 
                inplace=True)
    return df

def _add_ab_col(df: pd.DataFrame) -> pd.DataFrame:
    for rel in ['pro', 'pat', 'mat']:
        df[f'AB({rel})'] = df[f'AD({rel})'] / df[f'DP({rel})']
    return df

def _split_ggmacan_col(df: pd.DataFrame) -> pd.DataFrame:
    df = pd.concat(
        [df, df['GGM(AC/AN)'].str.split('/', expand=True)], axis=1)
    df = df.astype({0: float, 1: float})
    df.rename(columns={0: 'GGM(AC)', 1: 'GGM(AN)'}, inplace=True)

    return df


def _add_ggmaf_col(df: pd.DataFrame) -> pd.DataFrame:
    df['GGM(AF)'] = df['GGM(AC)'] / df['GGM(AN)']
    
    return df

def _rename_cols(df: pd.DataFrame) -> pd.DataFrame:
    rename_dict: dict = {
        'Chr': 'CHROM',
        'Position': 'POS',
        'Ref': 'REF',
        'Alt': 'ALT'
        }
    df.rename(columns=rename_dict, inplace=True)

    return df

def generate_variant_id(df: pd.DataFrame) -> pd.DataFrame:
    df['variant_id'] = df['Chr'] + ':' + \
                       df['Position'] + '-' + df['Ref'] + '-' + df['Alt']
    return df

def add_qc_filter(row, thresholds: Thresholds) -> pd.DataFrame:    
    if row['GQ(pro)'] < thresholds.GQ:
        return '.'
    else:
        if row['Vtype'] == 'homo':
            if row['AB(pro)'] >= (1 - thresholds.AB_hom):
                return 'PASS'
            else:
                return '.'
        else:
            if thresholds.AB_het <= row['AB(pro)'] <= (1 - thresholds.AB_het):
                return 'PASS'
            else:
                return '.'

def add_insilico_filter(row, thresholds: Thresholds) -> str:
    if ((row['SIFT'] >= thresholds.SIFT) 
        and (row['PolyPhen-2'] <= thresholds.PP2)
        and (row['CADD'] < thresholds.CADD)
        and (row['MutationTaster'] in thresholds.MutationTaster)):
        return 'FAIL' 
    else:
        return 'PASS'

def add_identified_filter(df: pd.DataFrame) -> pd.DataFrame:
    df.loc[
        df['Analysis status'] != 'Identified', 
        'Not_Identified_FILTER'] = 'PASS'

    return df

def add_ggmmaf_filter(df: pd.DataFrame, thresholds: Thresholds) -> pd.DataFrame:
    df.loc[
        df['GGM(AF)'] < thresholds.GGM_AF,
        'GGM_FILTER'] = 'PASS'
    
    return df
    
def add_hard_filter(df: pd.DataFrame) -> pd.DataFrame:
    df.loc[
        (
            (df['QC_FILTER'] == 'PASS')
            & (df['Not_Identified_FILTER'] == 'PASS')
            & (df['insilico_FILTER'] == 'PASS')
            & (df['GGM_FILTER'] == 'PASS')
        ),
        'HARD_FILTER'] = 'PASS'
    
    return df

# Pre-processing functions
def preprocess(df: pd.DataFrame) -> pd.DataFrame:
    df = _split_qc_col(df)
    df = _add_ab_col(df)
    df = _split_ggmacan_col(df)
    df = _add_ggmaf_col(df)
    df = generate_variant_id(df)
    df.replace({'SIFT': '-', 'PolyPhen-2': '-', 'CADD': '-'}, np.nan, inplace=True)
    df = df.astype(
        {'SIFT': float, 'PolyPhen-2': float, 'CADD': float, 
         'GQ(pro)': float, 'AB(pro)': float, 'GGM(AF)': float}
        )

    return df

def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    # Drop the 1st column
    df = df.drop(columns=df.columns[0])
    
    # Rename the columns
    df = _rename_cols(df)

    # Reorder the columns
    reodered_columns = [
        'HARD_FILTER', 'Gene', 'Transcript', 'Family', 'Sample', 'Disease', 
        'Vtype', 'variant_id', 'Amino acid change2', 'Effect',  
        'Distance', 'SIFT', 'PolyPhen-2', 'MutationTaster', 'CADD',
        'gnomAD(AF)', 'ExAC(AF)', 'ToMMo3.5K(AF)', 'GGM(AF)', 'JPNCTL(SC)', 
        'GGM(AC)', 'gnomAD(AC)', 'ToMMo3.5K(AC)',  
        'ID(pro)', 'AS(pro)', 'GT(pro)', 'GQ(pro)', 'DP(pro)', 'AD(pro)', 'AB(pro)',
        'ID(pat)', 'AS(pat)', 'GT(pat)', 'GQ(pat)', 'DP(pat)', 'AD(pat)', 'AB(pat)',
        'ID(mat)', 'AS(mat)', 'GT(mat)', 'GQ(mat)', 'DP(mat)', 'AD(mat)', 'AB(mat)',
        'Impact', 'QC_FILTER', 'Not_Identified_FILTER', 'insilico_FILTER',
        'Analysis status', 'Identified gene', 'Variant description', 
        'CHROM', 'POS', 'REF', 'ALT'
        ]
    
    df = df[reodered_columns]

    return df

def postprocess(df: pd.DataFrame) -> pd.DataFrame:
    df.replace(np.nan, '.', inplace=True)

    return df


def main(input_path: pathlib.Path):
    df = pd.read_csv(input_path, header=0, dtype=str)
    df = preprocess(df)
    df = add_identified_filter(df)
    df = add_ggmmaf_filter(df, Thresholds())
    df['QC_FILTER'] = df.apply(add_qc_filter, args=(Thresholds(),), axis=1)
    df['insilico_FILTER'] = df.apply(add_insilico_filter, args=(Thresholds(),), axis=1)
    df = add_hard_filter(df)
    df = reorder_columns(df)
    df = postprocess(df)

    df.to_excel(f"{input_path.parent}/{input_path.stem}.parsed.xlsx", index=False)


if __name__ == '__main__':
    input_path = pathlib.Path(sys.argv[1])
    main(input_path=input_path)