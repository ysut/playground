import pandas
import os

class Output:
    def __init__(self, outputDir, proband_id):
        self.outputDir = outputDir
        self.proband_id = proband_id
        self.saving_dir = f'{self.outputDir}/{self.proband_id}_results/'
        os.makedirs(self.saving_dir, exist_ok=True)
    
    def to_tsv(self, dataframe, fileName, index):
        dataframe.to_csv(f'{self.saving_dir}/{self.proband_id}_{fileName}.tsv', sep='\t', index=index)

    def to_excel(self, dataframe, fileName):
        dataframe.to_excel(f'{self.saving_dir}/{self.proband_id}_{fileName}.xlsx', index=False)
        return f'{self.saving_dir}/{self.proband_id}_{fileName}.xlsx'

    def to_stdout(self, countlist):
        df_count_summary = pandas.DataFrame(data=
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

