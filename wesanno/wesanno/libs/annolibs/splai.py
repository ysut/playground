from tqdm import tqdm
import glob
import pandas as pd
import os
import numpy as np
import pysam

class SplAI:
    def __init__(self, root_path):
        self._root_path = root_path
        self._snvFile = glob.glob(f'{self._root_path}/SpliceAI/*.snv.*.vcf.gz')
        self._indelFile = glob.glob(f'{self._root_path}/SpliceAI/*.indel.*.vcf.gz')
        
    def _pre_processing(self, dataframe):
        if 'variant_id' in dataframe.columns:
            return dataframe
        else:
            dataframe['varinat_id'] = dataframe['CHROM'] + '-' + dataframe['POS'] + '-' + dataframe['REF'] + '-' + dataframe['ALT']
            return dataframe
    
    def annotation(self, dataframe):
        if len(dataframe) != 0:
            dataframe = self._pre_processing(dataframe)
            dataframe = dataframe.assign(SpliceAI_INFO = 'NA|NA|NA|NA|NA|NA|NA|NA|NA|NA')
            self._chr_column = dataframe.columns.get_loc('CHROM')
            self._pos_column = dataframe.columns.get_loc('POS')
            self._splai_column = dataframe.columns.get_loc('SpliceAI_INFO')
            self._var_id_sr = dataframe['varinat_id']
            tbx_snv = pysam.TabixFile(self._snvFile[0])
            tbx_indel = pysam.TabixFile(self._indelFile[0])

            for index, variant_id in enumerate(tqdm(self._var_id_sr)):
                str_contig = dataframe.iat[index, self._chr_column]

                if str_contig == 'X' or str_contig == 'Y':
                    query_contig = str_contig
                    pass
                else:
                    query_contig = int(str_contig)
                    pass

                query_pos_start = int(dataframe.iat[index, self._pos_column]) -1
                query_pos_end = query_pos_start + 2

                for row in tbx_snv.fetch(query_contig, query_pos_start, query_pos_end, parser=pysam.asVCF()):
                    pos_snv = row.pos + 1
                    splai_id = f'{row.contig}-{pos_snv}-{row.ref}-{row.alt}'
                    if variant_id == splai_id :
                        dataframe.iat[index, self._splai_column] = row.info.replace('SpliceAI=', '')
                        break
                    else:
                        pass
                
                if dataframe.iat[index, self._splai_column]:
                    pass
                else:
                    for row in tbx_indel.fetch(query_contig, query_pos_start, query_pos_end, parser=pysam.asVCF()):
                        pos_indel = row.pos + 1
                        splai_id = f'{row.contig}-{pos_indel}-{row.ref}-{row.alt}'
                        if variant_id == splai_id :
                            dataframe.iat[index, self._splai_column] = row.info.replace('SpliceAI=', '')
                            break
                        else:
                            pass
            
            dataframe = pd.concat([dataframe, dataframe['SpliceAI_INFO'].str.split('|', expand=True)], axis=1)
            dataframe.rename(columns={0: 'splai_alt',
                                    1: 'splai_symbol',
                                    2: 'acp.gain',
                                    3: 'acp.loss',
                                    4: 'dnr.gain',
                                    5: 'dnr.loss',
                                    6: 'acp.gain_pos',
                                    7: 'acp.loss_pos',
                                    8: 'dnr.gain_pos',
                                    9: 'dnr.loss_pos'}, 
                                    inplace=True)
            dataframe.drop(columns=['SpliceAI_INFO', 'splai_alt', 'splai_symbol'], inplace=True)
            pass
        else:
            pass

        splai_scores_columns = ['acp.gain', 'acp.loss', 'dnr.gain', 'dnr.loss']
        dataframe = dataframe.dropna(subset=splai_scores_columns)
        dataframe['max_splai'] = dataframe[splai_scores_columns].max(axis=1)
        dataframe['max_splai'] = dataframe['max_splai'].replace('NA', np.nan)


        return dataframe





