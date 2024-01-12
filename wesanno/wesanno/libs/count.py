from .filtering import Filter
import pandas as pd

class FilterCount(Filter):
    def __init__(self, dataframe, config_maf, *args):
        
        self._countAD = []
        self._countHm = []
        self._countCH = []
        self._countXL = []
        self.countlist = [self._countAD, self._countHm, self._countCH, self._countXL]
        self.dataframe = dataframe
        self.exac_cutoff = config_maf['exac']
        self.esp_cutoff = config_maf['esp']
        self.hgvd_cutoff = config_maf['hgvd']
        self.tommo_cutoff = config_maf['tommo']
        self.inhouse_cutoff = config_maf['inhouse']
        self.df = super().pre_genotyping(self.dataframe, *args)
        self.df = super().pre_processing(self.df)
        self._dfA, self._dfS = super().chrom_filter(self.df)
        self._dfA, self._dfS = super().pre_inhouse(self._dfA, self._dfS)
        super().count(self._dfA, self._dfA, self._dfA, self._dfS)

    # def all_filter(self, ):
    #     df = super().pre_processing(df)
    
    #     fltr = Filter(config['Minor_allele_frequency'])
    #     df = fltr.pre_genotyping(df, *(pickup.pickup_dnaID()))
    #     df = fltr.pre_processing(df)
    #     dfA, dfS = fltr.chrom_filter(df)
    #     dfA, dfS = fltr.pre_inhouse(dfA, dfS)
    #     dfA, dfS = fltr.cast(dfA, dfS)
    #     fltr.count(dfA, dfA, dfA, dfS, dfS)

    #     ## Filtering inhouse
    #     df_AD = fltr.inhouse_filter(dfA, 'ad')
    #     df_Hm = fltr.inhouse_filter(dfA, 'homo')
    #     df_CH = fltr.inhouse_filter(dfA, 'ch')
    #     df_XL = fltr.inhouse_filter(dfS, 'xlink')
    #     fltr.count(df_AD, df_Hm, df_CH, df_XL)

    #     ## Filtering MAF
    #     df_AD = fltr.maf_filter(df_AD, 'ad')
    #     df_Hm = fltr.maf_filter(df_Hm, 'homo')
    #     df_CH = fltr.maf_filter(df_CH, 'ch')
    #     df_XL = fltr.maf_filter(df_XL, 'xlink')
    #     fltr.count(df_AD, df_Hm, df_CH, df_XL)

    #     ## Excluding HLA/MUC
    #     df_AD = fltr.excluding_hlamuc(df_AD)
    #     df_Hm = fltr.excluding_hlamuc(df_Hm)
    #     df_CH = fltr.excluding_hlamuc(df_CH)
    #     df_XL = fltr.excluding_hlamuc(df_XL)
    #     fltr.count(df_AD, df_Hm, df_CH, df_XL)

    #     ## Excluding Synoynmous
    #     df_AD = fltr.excluding_synonymous(df_AD)
    #     df_Hm = fltr.excluding_synonymous(df_Hm)
    #     df_CH = fltr.excluding_synonymous(df_CH)
    #     df_XL = fltr.excluding_synonymous(df_XL)
    #     fltr.count(df_AD, df_Hm, df_CH, df_XL)

    #     ## GT filter
    #     dataframe_dic = {'df_AD': df_AD, 'df_Hm': df_Hm, 'df_CH': df_CH, 'df_XL': df_XL}
    #     dataframe_tuple = fltr.gt_filter(args['mode'], **dataframe_dic)
    #     df_AD, df_Hm, df_CH, df_XL = dataframe_tuple
    #     countlist = fltr.count(df_AD, df_Hm, df_CH, df_XL)

    #     ### Stdoutput Summary file
    #     out.to_stdout(countlist)
