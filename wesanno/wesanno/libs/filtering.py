import pandas as pd

class Filter:
    def __init__(self, dataframe, mode, config_maf, *args):
        self.df = dataframe
        self.mode = mode
        self.countAD = []
        self.countHm = []
        self.countCH = []
        self.countXL = []
        self.countlist = [self.countAD, self.countHm, self.countCH, self.countXL]
        self.exac_cutoff = config_maf['exac']
        self.esp_cutoff = config_maf['esp']
        self.hgvd_cutoff = config_maf['hgvd']
        self.tommo_cutoff = config_maf['tommo']
        self.inhouse_cutoff = config_maf['inhouse']
        self.sample_ids = args # self.sample_ids is a tuple

    def count(self, *args):
        for model, dataframe in zip(self.countlist, args):
            model.append(len(dataframe))
        
    def pre_processing(self, dataframe):
        dataframe['HGVD_AF'] = dataframe['HGVD'].str.extract('(\.|[0-1]\.\d{,6})')
        dataframe['ToMMo3.5KJPN_AF'] = dataframe['snp20171005_tommo3.5k_passed'].str.extract('(\.|[0-1]\.\d{0,})')
        dataframe['ExAC_ALL'] = dataframe['ExAC_ALL'].replace('.','0')
        dataframe['HGVD_AF'] = dataframe['HGVD_AF'].replace('.','0')
        dataframe['esp6500siv2_all'] = dataframe['esp6500siv2_all'].replace('.','0')
        dataframe['ToMMo3.5KJPN_AF'] = dataframe['ToMMo3.5KJPN_AF'].replace('.','0')
        return dataframe
    
    def pre_genotyping(self, dataframe, *args):
        sample_num = len(args)
        if sample_num == 1:
            proband_id = args[0]
            dataframe['GT_Pro'] = dataframe[proband_id].str.extract('(\./\.|\d/\d)')
            pass
        elif sample_num == 3:
            proband_id, father_id, mother_id = args
            dataframe['GT_Pro'] = dataframe[proband_id].str.extract('(\./\.|\d/\d)')
            dataframe['GT_Fa'] = dataframe[father_id].str.extract('(\./\.|\d/\d)')
            dataframe['GT_Mo'] = dataframe[mother_id].str.extract('(\./\.|\d/\d)')
            pass
        elif sample_num == 4:
            proband_id, father_id, mother_id, sibling_id = args
            dataframe['GT_Pro'] = dataframe[proband_id].str.extract('(\./\.|\d/\d)')
            dataframe['GT_Fa'] = dataframe[father_id].str.extract('(\./\.|\d/\d)')
            dataframe['GT_Mo'] = dataframe[mother_id].str.extract('(\./\.|\d/\d)') 
            dataframe['GT_Sib'] = dataframe[sibling_id].str.extract('(\./\.|\d/\d)') 
            pass
        elif sample_num == 2:
            proband_id, parent_id = args
            dataframe['GT_Pro'] = dataframe[proband_id].str.extract('(\./\.|\d/\d)')
            dataframe['GT_Par'] = dataframe[parent_id].str.extract('(\./\.|\d/\d)')
            pass
        else:
            print('Erorr')
            pass
        return dataframe

    def chrom_filter(self, dataframe):
        dfA = dataframe[(dataframe['Chr'] != 'X') & (dataframe['Chr'] != 'Y')].copy()
        dfS = dataframe[(dataframe['Chr'] == 'X') | (dataframe['Chr'] == 'Y')].copy()
        return dfA, dfS
    
    def pre_inhouse(self, dataframeA, dataframeS):
        dataframeA['Inhouse'] = dataframeA['INFO'].str.extract('JpnMutation=(p|\+|\-)')
        dataframeA['Inhouse_575'] = dataframeA['INFO'].str.extract(':575:(\d{,3})')
        dataframeS['Inhouse'] = dataframeS['INFO'].str.extract('JpnMutation=(p|\+|\-)')
        dataframeS['Inhouse_M'] = dataframeS['INFO'].str.extract('281,294:(\d{,3})')
        dataframeS['Inhouse_F'] = dataframeS['INFO'].str.extract('281,294:\d{,3},(\d{,3})')
        return dataframeA, dataframeS

    def cast(self, dataframeA, dataframeS):
        dataframeA = dataframeA.astype({'Inhouse_575': 'int', 'ExAC_ALL': 'float', 
                                        'HGVD_AF': 'float', 'esp6500siv2_all': 'float', 
                                        'ToMMo3.5KJPN_AF': 'float'})
        dataframeS = dataframeS.astype({'Inhouse_M': 'int','Inhouse_F': 'int', 
                                        'ExAC_ALL': 'float','HGVD_AF': 'float', 
                                        'esp6500siv2_all': 'float', 'ToMMo3.5KJPN_AF': 'float'})
        return dataframeA, dataframeS

    def inhouse_filter(self, dataframe, model):
        chr_column = dataframe.columns.get_loc('CHROM')
        if dataframe.iat[0, chr_column] == 'X' or dataframe.iat[0, chr_column] == 'Y':
            dataframe = dataframe[(dataframe['Inhouse_M'] <= self.inhouse_cutoff[model]) & 
                                  (dataframe['Inhouse_F'] <= self.inhouse_cutoff[model])]
            pass
        else:
            dataframe = dataframe[dataframe['Inhouse_575'] <= self.inhouse_cutoff[model]]
            pass
        return dataframe

    def maf_filter(self, dataframe, model):
        dataframe = dataframe[dataframe['ExAC_ALL'] < self.exac_cutoff[model]]
        dataframe = dataframe[dataframe['esp6500siv2_all'] < self.esp_cutoff[model]]
        dataframe = dataframe[dataframe['HGVD_AF'] < self.hgvd_cutoff[model]]
        dataframe = dataframe[dataframe['ToMMo3.5KJPN_AF'] < self.tommo_cutoff[model]]
        return dataframe

    def excluding_hlamuc(self, dataframe):
        dataframe = dataframe[~dataframe['Gene.refGene'].str.contains('HLA|MUC')]
        return dataframe

    def excluding_synonymous(self, dataframe):
        df = self.pre_genotyping(self.df, *self.sample_ids)
        df = self.pre_processing(df)
        dataframe = dataframe[~((dataframe['ExonicFunc.refGene'] == 'synonymous SNV') &
                               (dataframe['Func.refGene'] == 'exonic'))]
        return dataframe
    
    def extract_synonymous(self, dataframe):
        dataframe = dataframe[dataframe['ExonicFunc.refGene'] == 'synonymous SNV']
        return dataframe

    def _hetero_count(self, dataframe):
        series = dataframe.groupby('Gene.refGene')['Gene.refGene'].transform('count').rename('Hetero_count')
        dataframe = pd.concat([dataframe, series], axis=1)
        return dataframe

    def _gt_xl_filter(self, dataframe):
        dataframe = dataframe[~dataframe['GT_Pro'].str.contains('0/0|\./\.')]
        return dataframe
    
    def _gt_ch_filter(self, dataframe):
        dataframe = self._hetero_count(dataframe)
        dataframe = dataframe[dataframe['GT_Pro'] == '0/1']
        dataframe = dataframe[dataframe['Hetero_count'] >= 2]
        return dataframe

    def _gt_hm_filter(self, dataframe, mode):
        if mode == 'proband' or mode == 'duo':
            dataframe = dataframe[(dataframe['GT_Pro'] == '1/1')]
            return dataframe
        else:
            dataframe = dataframe[(dataframe['GT_Pro'] == '1/1') & 
                                  (dataframe['GT_Fa'] != '1/1') &
                                  (dataframe['GT_Mo'] != '1/1')]
            return dataframe

    def _gt_ad_filter(self, dataframe, mode):
        if mode == 'proband':
            dataframe = self._hetero_count(dataframe)
            dataframe = dataframe[dataframe['GT_Pro'] == '0/1']
            dataframe = dataframe[dataframe['Hetero_count'] == 1]
            return dataframe
        elif mode == 'duo':
            dataframe = self._hetero_count(dataframe)
            dataframe = dataframe[(~dataframe['GT_Pro'].str.contains('0/0|\./\.')) &
                                  (dataframe['GT_Par'].str.contains('0/0|\./\.'))]
            dataframe = dataframe[dataframe['Hetero_count'] == 1]
            return dataframe
        else:
            dataframe = dataframe[(~dataframe['GT_Pro'].str.contains('0/0|\./\.')) &
                                  (dataframe['GT_Fa'].str.contains('0/0|\./\.')) &
                                  (dataframe['GT_Mo'].str.contains('0/0|\./\.'))]
            return dataframe

    def gt_filter(self, mode, **kwargs):
        df_AD = self._gt_ad_filter(kwargs['df_AD'], mode)
        df_Hm = self._gt_hm_filter(kwargs['df_Hm'], mode)
        df_CH = self._gt_ch_filter(kwargs['df_CH'])
        df_XL = self._gt_xl_filter(kwargs['df_XL'])
        return df_AD, df_Hm, df_CH, df_XL

    def gt_snvcnv_filter(self, dataframe):
        dataframe = dataframe[(dataframe['GT_Pro'] == '0/1')]
        return dataframe


    def all_filter(self):
        # print(*self.sample_ids)
        # print(type(self.sample_ids))
        
        # idCount = len(self.sample_ids)
        # print(f"Number of sample: {idCount}")
        # for i in range(idCount):
        #     print(self.sample_ids[i])

        self.df = self.pre_genotyping(self.df, *self.sample_ids)
        self.df = self.pre_processing(self.df)
        dfA, dfS = self.chrom_filter(self.df)
        dfA, dfS = self.pre_inhouse(dfA, dfS)
        dfA, dfS = self.cast(dfA, dfS)
        self.count(dfA, dfA, dfA, dfS)

        df_AD = self.inhouse_filter(dfA, 'ad')
        df_Hm = self.inhouse_filter(dfA, 'homo')
        df_CH = self.inhouse_filter(dfA, 'ch')
        df_XL = self.inhouse_filter(dfS, 'xlink')
        self.count(df_AD, df_Hm, df_CH, df_XL)

        ## Filtering MAF
        df_AD = self.maf_filter(df_AD, 'ad')
        df_Hm = self.maf_filter(df_Hm, 'homo')
        df_CH = self.maf_filter(df_CH, 'ch')
        df_XL = self.maf_filter(df_XL, 'xlink')
        self.count(df_AD, df_Hm, df_CH, df_XL)

        ## Excluding HLA/MUC
        df_AD = self.excluding_hlamuc(df_AD)
        df_Hm = self.excluding_hlamuc(df_Hm)
        df_CH = self.excluding_hlamuc(df_CH)
        df_XL = self.excluding_hlamuc(df_XL)
        self.count(df_AD, df_Hm, df_CH, df_XL)

        df_AD = self.excluding_synonymous(df_AD)
        df_Hm = self.excluding_synonymous(df_Hm)
        df_CH = self.excluding_synonymous(df_CH)
        df_XL = self.excluding_synonymous(df_XL)
        self.count(df_AD, df_Hm, df_CH, df_XL)
        df_SC = df_Hm

        ## GT filter
        dataframe_dic = {'df_AD': df_AD, 'df_Hm': df_Hm, 
                         'df_CH': df_CH, 'df_XL': df_XL}
        dataframe_tuple = self.gt_filter(self.mode, **dataframe_dic)
        df_AD, df_Hm, df_CH, df_XL = dataframe_tuple
        self.count(df_AD, df_Hm, df_CH, df_XL)

        return df_AD, df_Hm, df_CH, df_XL, df_SC
