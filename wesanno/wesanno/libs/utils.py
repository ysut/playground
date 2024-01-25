import os
import pathlib2
import pandas as pd
import toml

from dataclasses import dataclass


def load_config(path_to_tomlFile) -> dict:
    with open(path_to_tomlFile) as f:
        config: dict = toml.load(f)
    
    return config

def dfs_to_excel(dfs: dataclass, output_xlsx: str) -> None:
    sheet_names = ['AD', 'Homo', 'CH', 'XL']
    with pd.ExcelWriter(output_xlsx) as writer:
        for df, sheet_name in zip([dfs.AD, dfs.Hm, dfs.CH, dfs.XL], sheet_names):
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    
    return None


class OutputSettings:
    def __init__(self, args, mode_samples_info):
        self.proband_id = mode_samples_info.proband_id
        self.mode = mode_samples_info.mode
        if args['output']:
            self.output_dir = args['output']
        else:
            path_input_file = pathlib2.Path(args['input'])
            self.output_dir = path_input_file.parent
        
        self.saving_dir = f'{self.output_dir}/{self.proband_id}-{self.mode}_results'
        os.makedirs(self.saving_dir, exist_ok=True)
    
    def get_saving_file_path(self) -> str:
        return f'{self.saving_dir}/{self.proband_id}-{self.mode}.tsv'


    # 以下は使わないかもしれん

    def to_tsv(self, dataframe, fileName, index):
        dataframe.to_csv(f'{self.saving_dir}/{self.proband_id}_{fileName}.tsv', sep='\t', index=index)

    def to_excel(self, dataframe, fileName):
        dataframe.to_excel(f'{self.saving_dir}/{self.proband_id}_{fileName}.xlsx', index=False)
        return f'{self.saving_dir}/{self.proband_id}_{fileName}.xlsx'

    def to_stdout(self, countlist):
        df_count_summary = pd.DataFrame(data=
                                {'AD model': countlist[0],
                                    'Homo model': countlist[1],
                                    'CH model': countlist[2],
                                    'X-linked model': countlist[3]},
                                    index=
                                    ['Chromosome',
                                    'In-house filter', 
                                    'MAF filter', 
                                    'Excluding HLA/MUC', 
                                    'Excluding exonic synonymous',
                                    'GT filter']
                                    )
        self.to_tsv(df_count_summary, 'CountSummary', True)
        print(df_count_summary)
        print()
        return df_count_summary

