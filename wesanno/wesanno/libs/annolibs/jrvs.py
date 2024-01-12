from tqdm import tqdm
import glob
import pandas as pd
import pysam

class Jrvs:

    def __init__(self, root_path):
        self._root_path = root_path
        self._snvFile = glob.glob(f'{self._root_path}/SpliceAI/*.snv.*.vcf.gz')
        self._indelFile = glob.glob(f'{self._root_path}/SpliceAI/*.indel.*.vcf.gz')

    def annotation(self, dataframe):
        dataframe['JARVIS'] = None
        dataframe['%JARVIS'] = None
        self._chr_column = dataframe.columns.get_loc('CHROM')
        self._pos_column = dataframe.columns.get_loc('POS')
        jrvs_column = dataframe.columns.get_loc('JARVIS')
        pjrvs_column = dataframe.columns.get_loc('%JARVIS')

        for index in tqdm(range(len(dataframe))): 
            query_contig = dataframe.iat[index, self._chr_column]

            if query_contig == 'X' or query_contig == 'Y':
                break
            else:
                pass

            path_to_jrvsFile = glob.glob(f'{self._root_path}/JARVIS/*.chr{query_contig}.*.gz')                       
            tbx_jrvs = pysam.TabixFile(path_to_jrvsFile[0])
            query_start = int(dataframe.iat[index, self._pos_column])
            query_end = query_start + 1
            query_id = f'chr{query_contig}-{query_start}'

            for row in tbx_jrvs.fetch(f'chr{query_contig}', query_start, query_end, parser=pysam.asBed()):
                jrvs_id = f'{row.contig}-{row.start}'

                if query_id == jrvs_id:
                    dataframe.iat[index, jrvs_column] = row[3]
                    dataframe.iat[index, pjrvs_column] = row[4]
                    break
                else:
                    pass

        return dataframe

