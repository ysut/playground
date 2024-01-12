import pandas as pd

class DFprocess:
    def __init__(self, non_xhmm_dfs: list, xhmm_df=None):
        self.dfs = non_xhmm_dfs
        self.df_xhmm = xhmm_df

    def _replace_pred(self):
        for df in self.dfs:
            # Binary
            df.replace({'MetaSVM_pred': {'D': 4, 'T': 2}}, inplace=True)
            df.replace({'MetaLR_pred': {'D': 4, 'T': 2}}, inplace=True)
            df.replace({'fathmm-MKL_coding_pred': {'D': 4, 'N': 2}}, inplace=True)
            df.replace({'FATHMM_pred': {'D': 4, 'T': 2}}, inplace=True)
            # Tertiary
            df.replace({'LRT_pred': {'D': 4, 'U': 3, 'N': 2}}, inplace=True)
            # Quaternary
            df.replace({'MutationAssessor_pred': {'H': 4, 'M': 3, 'L': 2, 'N': 2}}, inplace=True)

    def _cast_for_non_xhmm(self):
        casted_dfs = []

        full_cast_columns = ['esp6500siv2_all', 'SIFT_score', 'Polyphen2_HVAR_score', 
                             'MutationTaster_score', 'CADD_phred', 'DANN_score', 
                             'MutationAssessor_score', 'FATHMM_score', 'PROVEAN_score',
                             'GERP++_RS', 'phyloP7way_vertebrate', 
                             'phyloP20way_mammalian', 'phastCons7way_vertebrate',
                             'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'gerp++gt2',
                             'pLI', 'pRec', 'LOEUF', 'syn_z', 'mis_z', 
                             'PROVEAN_score', 'HGVD_AF', 'ToMMo3.5KJPN_AF', 
                             'Hetero_count', 'Start', 'End', 'Inhouse_575', 'MetaSVM_pred', 
                             'MetaLR_pred', 'fathmm-MKL_coding_pred', 'LRT_pred', 
                             'MutationAssessor_pred', 'Start', 'End', 'max_splai'] 
        
        for df in self.dfs:
            preproc_float_columns = []
            for cast_colmun in full_cast_columns:
                if cast_colmun in df.columns:
                    preproc_float_columns.append(cast_colmun)
                else:
                    pass
        
            # make dictionary for casting
            astype_float_dict = {column: float for column in preproc_float_columns}
            for column in preproc_float_columns:
                df[column] = df[column].mask(df[column] == '.') 
        
            # cast and append to list
            df = df.astype(astype_float_dict)
            casted_dfs.append(df)
        
        return casted_dfs

    def _cast_for_xhmm(self):
        # preproc_float_columns = ['KB', 'CHR', 'MID_BP',
        #                          'NUM_TARG', 'Q_EXACT', 'Q_SOME',
        #                          'Q_NON_DIPLOID', 'Q_START', 'Q_STOP',
        #                          'MEAN_RD', 'MEAN_ORIG_RD', 'SegDup',
        #                          'StrVar', 'Haplo_Insuff_Trip', 
        #                          'Decipher', 'Omimgene']
        preproc_float_columns = ['KB', 'MID_BP',
                                 'NUM_TARG', 'Q_EXACT', 'Q_SOME',
                                 'Q_NON_DIPLOID', 'Q_START', 'Q_STOP',
                                 'MEAN_RD', 'MEAN_ORIG_RD', 'SegDup',
                                 'StrVar', 'Haplo_Insuff_Trip', 
                                 'Decipher']
        astype_float_dict = {column: float for column in preproc_float_columns}
        for column in preproc_float_columns:
            self.df_xhmm[column] = self.df_xhmm[column].mask(self.df_xhmm[column] == '.') 
        df_xhmm = self.df_xhmm.astype(astype_float_dict)

        return df_xhmm

    def preprocess_for_non_xhmm(self):
        self._replace_pred()
        casted_dfs = self._cast_for_non_xhmm()
        return casted_dfs   # retrun a list of processed dataframes

    def preprocess_for_xhmm(self):
        casted_df = self._cast_for_xhmm()
        return casted_df   # retrun a dataframe

    def rearrange_column_order(self):
        pass