import pandas as pd

'''
input: df

1. annotate evaluation summary to dfs (df_AD, df_Hm, df_CH, df_XL) 
    1.1. prediction tools
         SIFT, PP2_hvar, MutationTaster, CADD
    1.2.  
3. setting column regions for openpyxl
4.

'''


class Evaluation:
    def __init__(self, dataframes: list):
        self.dfs = dataframes
        pass

    def eval_inheritance_duo_ProPar(self):
        for df in self.dfs:
            df.loc[(df['GT_Pro'] == '0/1')
                    & ((df['GT_Par'] == '0/0') | (df['GT_Par'] == './.')),
                    'Inheritance'] = 'DeNovo?'
            df.loc[((df['GT_Pro'] == '0/1') | (df['GT_Pro'] == '1/1'))
                    & ((df['GT_Par'] == '0/1') | (df['GT_Par'] == '1/1') | (df['GT_Par'] == '1/2')),
                    'Inheritance'] = 'inh_Par?'

    def eval_inheritance_trio(self):
        for df in self.dfs:
            df.loc[(df['GT_Pro'] == '0/1')
                    & ((df['GT_Fa'] == '0/0') | (df['GT_Fa'] == './.'))
                    & ((df['GT_Mo'] == '0/0') | (df['GT_Fa'] == './.')),
                    'Inheritance'] = 'DeNovo'
            df.loc[(df['GT_Pro'] == '0/1')
                    & (df['GT_Fa'] == '0/1')
                    & (df['GT_Mo'] == '0/0'),
                    'Inheritance'] = 'inh_Fa'
            df.loc[(df['GT_Pro'] == '0/1')
                    & (df['GT_Fa'] == '0/1')
                    & (df['GT_Mo'] == './.'),
                    'Inheritance'] = 'inh_Fa?'
            df.loc[(df['GT_Pro'] == '0/1')
                    & (df['GT_Fa'] == '0/0')
                    & (df['GT_Mo'] == '0/1'),
                    'Inheritance'] = 'inh_Mo'
            df.loc[(df['GT_Pro'] == '0/1')
                    & (df['GT_Fa'] == './.')
                    & (df['GT_Mo'] == '0/1'),
                    'Inheritance'] = 'inh_Mo?'
            df.loc[(df['GT_Pro'] == '0/1')
                    & (df['GT_Fa'] == '1/1')
                    & (df['GT_Mo'] == '1/1'),
                    'Inheritance'] = 'inh_Fa?Mo?'
        
    def sibling_analysis(self):
        for df in self.dfs:
            df.loc[(df['GT_Pro'] == '0/1') 
                   & (df['GT_Sib'] == '0/1'),
                   'QuadAnalysis'] = 'Pro+/Sib+'
            df.loc[(df['GT_Pro'] == '0/1') 
                   & (df['GT_Sib'] == '0/0') ,
                   'QuadAnalysis'] = 'Pro+/Sib-'
            df.loc[(df['GT_Pro'] == '0/1') 
                   & (df['GT_Sib'] == './.') ,
                   'QuadAnalysis'] = 'Pro+/Sib?'
    
    def eval_clinsig_info(self):
        for df in self.dfs:
            df.loc[(df['CLINSIG'] == '0/1') 
                   & (df['GT_Sib'] == '0/1'),
                   'QuadAnalysis'] = 'Pro+/Sib+'
            df.loc[(df['GT_Pro'] == '0/1') 
                   & (df['GT_Sib'] == '0/0') ,
                   'QuadAnalysis'] = 'Pro+/Sib-'
            df.loc[(df['GT_Pro'] == '0/1') 
                   & (df['GT_Sib'] == './.') ,
                   'QuadAnalysis'] = 'Pro+/Sib?'
        
    def add_comment_if_string_contains(self, df, column_name, string_to_search, comment_column_name, comment):
        df[comment_column_name] = df.apply(lambda row: comment if string_to_search in row[column_name] else row[column_name], axis=1)


    def _missense_eval(self):
        for df in self.dfs:
            df.loc[(df['MutationTaster_pred'] == 'D') | 
                        (df['MutationTaster_pred'] == 'A'), 
                        'MutationTasterSummary'] = 'DorA'
            df.loc[(df['MutationTaster_pred'] == 'N') | 
                        (df['MutationTaster_pred'] == 'P'), 
                        'MutationTasterSummary'] = 'NorP'
            df.loc[(df['SIFT_score'] < 0.05) &
                        (df['Polyphen2_HVAR_score'] > 0.909) & 
                        (df['MutationTasterSummary'] == 'DorA'),  
                        'PredictionToolsSummary'] = 'Deleterious'
            df.loc[(df['SIFT_score'] < 0.05) &
                        (df['Polyphen2_HVAR_score'] > 0.909) & 
                        (df['MutationTasterSummary'] == 'DorA') &
                        (df['CADD_phred'] > 25),
                        'PredictionToolsSummary'] = 'Deleterious!'
            df.loc[(df['SIFT_score'] < 0.05) &
                        (df['Polyphen2_HVAR_score'] <= 0.909) &
                        (df['Polyphen2_HVAR_score'] > 0.447) & 
                        (df['MutationTasterSummary'] == 'DorA'),  
                        'PredictionToolsSummary'] = 'PossiblyDeleterious'
            df.loc[(df['SIFT_score'] >= 0.05) &
                        (df['Polyphen2_HVAR_score'] <= 0.447) & 
                        (df['MutationTasterSummary'] == 'NorP'),  
                        'PredictionToolsSummary'] = 'PossiblyBenign'
            

    def eval_clinsig():
        pass

    def eval_():
        pass


    def evaluatate_all(self):
        self._eval_missense()
        pass

class PrimaryEvaluation:
    def __init__(self, dataframe):
        self.df = dataframe
        pass


