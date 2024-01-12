import pandas as pd
import re
import numpy as np
import gzip

class Xhmm:
    def __init__(self, std, xx, xy):
        self.std = pd.read_table(std, header=0, dtype=str)
        self.xx = pd.read_table(xx, header=0, dtype=str)
        self.xy = pd.read_table(xy, header=0, dtype=str)
        self.std = self.std[(self.std['CHR'] != 'X') & (self.std['CHR'] != 'Y')]
        self.xx = self.xx[self.xx['CHR'] == 'X']
        self.xy = self.xy[self.xy['CHR'] == 'Y']

    def make_score_list(self, path_to_db):
        self.cnv_score_list = []
        with gzip.open(path_to_db, 'rt') as f:
            for buff in f:
                self.cnv_score_list.append(buff.split())
        self.cnv_gene_set = set([row[0] for row in self.cnv_score_list])
        return self.cnv_score_list

    def _extract_samples(self, dataframe, *args):
        sample_num = len(args)
        if sample_num == 1:
            proband_id = args[0]
            dataframe = dataframe[dataframe['SAMPLE'] == proband_id]
        elif sample_num == 2:
            proband_id, parent_id = args
            dataframe = dataframe[(dataframe['SAMPLE'] == proband_id) |
                                  (dataframe['SAMPLE'] == parent_id)]
        elif sample_num == 3:
            proband_id, father_id, mother_id = args
            dataframe = dataframe[(dataframe['SAMPLE'] == proband_id) |
                                  (dataframe['SAMPLE'] == father_id) |
                                  (dataframe['SAMPLE'] == mother_id)]
        elif sample_num == 4:
            proband_id, father_id, mother_id, sibling_id = args
            dataframe = dataframe[(dataframe['SAMPLE'] == proband_id) |
                                  (dataframe['SAMPLE'] == father_id) |
                                  (dataframe['SAMPLE'] == mother_id) |
                                  (dataframe['SAMPLE'] == sibling_id)]
        else:
            print('IDs Error')
        return dataframe
    
    def __cast(self, dataframe):
        astype_dict = {'Q_SOME': 'int8', 'SegDup': 'float'}
        dataframe = dataframe.astype(astype_dict)
        return dataframe

    def _filter(self, dataframe, config):
        dataframe = self.__cast(dataframe)
        dataframe = dataframe[dataframe['Q_SOME'] >= int(config['XHMM']['qsome'])]
        dataframe = dataframe[dataframe['SegDup'] <= float(config['XHMM']['segdup'])]
        return dataframe  

    def _add_columns(self, dataframe):
        dataframe.insert(len(dataframe.columns), 'pHaplo', np.nan)
        dataframe.insert(len(dataframe.columns), 'pTriplo', np.nan)
        dataframe.insert(len(dataframe.columns), 'DosageSensitive', np.nan)
        return dataframe

    def annotate_and_evaluate_cnvs(self, dataframe, config, *args):
        dataframe = self._extract_samples(dataframe, *args)
        dataframe = self._filter(dataframe, config)

        if len(dataframe.index) == 0:
            return dataframe
        
        else:
            dataframe = self._add_columns(dataframe)
            cnv_col = dataframe.columns.get_loc('CNV')
            gene_col = dataframe.columns.get_loc('gene')
            phaplo_col = dataframe.columns.get_loc('pHaplo')
            ptriplo_col = dataframe.columns.get_loc('pTriplo')
            summary_col = dataframe.columns.get_loc('DosageSensitive')
            for index in range(len(dataframe)):
                phaplo_list = []
                ptriplo_list = []
                value = dataframe.iat[index,gene_col]
                genes = value.split(',')
                for gene in genes:
                    if gene not in self.cnv_gene_set:
                        phaplo_list.append('-')
                        ptriplo_list.append('-')
                        pass
                    else:
                        for data in self.cnv_score_list:
                            if gene == data[0]:
                                phaplo_list.append(data[1])
                                ptriplo_list.append(data[2])
                                if dataframe.iat[index,cnv_col] == 'DEL':
                                    if float(data[1]) >= config['Dosage_sensitivity']['phaplo']['precision']:
                                        dataframe.iat[index,summary_col] = 'Including Sensitive Gene'
                                    elif float(data[1]) >= config['Dosage_sensitivity']['phaplo']['recall']:
                                        dataframe.iat[index,summary_col] = 'Including Possibly Sensitive Gene'
                                elif dataframe.iat[index,cnv_col] == 'DUP':
                                    if float(data[2]) >= config['Dosage_sensitivity']['ptriplo']['precision']:
                                        dataframe.iat[index,summary_col] = 'Including Sensitive Gene' 
                                    elif float(data[2]) >= config['Dosage_sensitivity']['ptriplo']['recall']:
                                        dataframe.iat[index,summary_col] = 'Including Possibly Sensitive Gene' 
                                break
                            else:
                                pass

                if dataframe.iat[index,cnv_col] == 'DEL':
                    dataframe.iat[index,phaplo_col] = str(phaplo_list)      # Lists cannot be input
                elif dataframe.iat[index,cnv_col] == 'DUP':
                    dataframe.iat[index,ptriplo_col] = str(ptriplo_list)    # Lists cannot be input
                else:
                    pass   
            return dataframe
    
    def annotate_dosage_scores(self, dataframe):
        dataframe.insert(len(dataframe.columns), 'pHaplo', np.nan)
        dataframe.insert(len(dataframe.columns), 'pTriplo', np.nan)
        symbol_col = dataframe.columns.get_loc('Gene.refGene')
        phaplo_col = dataframe.columns.get_loc('pHaplo')
        ptriplo_col = dataframe.columns.get_loc('pTriplo')
        for i in range(len(dataframe)):
            value = dataframe.iat[i,symbol_col]
            genes = value.split(',')
            phaplo_list = []
            ptriplo_list = []
            for gene in genes:
                if gene not in self.cnv_gene_set:
                    phaplo_list.append('-')
                    ptriplo_list.append('-')
                    pass
                else:
                    for data in self.cnv_score_list:
                        if gene == data[0]:
                            phaplo_list.append(data[1])
                            ptriplo_list.append(data[2])
                            break
                        else:
                            pass
                dataframe.iat[i, phaplo_col] = phaplo_list
                dataframe.iat[i, ptriplo_col] = ptriplo_list
        return dataframe

    def __make_cnv_gene_set(self, series):
        gene_list = []
        for row in series:
            genes = row.split(',')
            gene_list.extend(genes)
        gene_set = set(gene_list)
        try:
            gene_set.remove('-') 
        except KeyError:
            pass
        return gene_set
    
    def _make_cnv_gene_sets(self, dataframe):
        dataframe_del = dataframe[dataframe['CNV'] == 'DEL']
        dataframe_dup = dataframe[dataframe['CNV'] == 'DUP']
        series_del = dataframe_del['gene']
        series_dup = dataframe_dup['gene']
        del_set = self.__make_cnv_gene_set(series_del)
        dup_set = self.__make_cnv_gene_set(series_dup)
        return del_set, dup_set
    
    def _union_set(self, dataframe_std, dataframe_xx, dataframe_xy):
        std_del_set, std_dup_set = self._make_cnv_gene_sets(dataframe_std)
        xx_del_set, xx_dup_set = self._make_cnv_gene_sets(dataframe_xx)
        xy_del_set, xy_dup_set = self._make_cnv_gene_sets(dataframe_xy)
        return std_del_set | xx_del_set | xy_del_set, std_dup_set | xx_dup_set | xy_dup_set

    def fetch_xhmm_genes(self, dataframe):
        del_set, dup_set = self._union_set(self.std, self.xx, self.xy)
        dataframe.insert(len(dataframe.columns), 'XHMM.info', np.nan)
        xhmminfo_col = dataframe.columns.get_loc('XHMM.info')
        gene_col = dataframe.columns.get_loc('Gene.refGene')
        for index in range(len(dataframe)):
            value = dataframe.iat[index,gene_col]
            genes = value.split(',')
            for gene in genes:
                if gene in del_set:
                    dataframe.iat[index,xhmminfo_col] = 'DEL gene'
                    break
                elif gene in dup_set:
                    dataframe.iat[index,xhmminfo_col] = 'Dup gene'
                    break
                else:
                    pass
                
        return dataframe
