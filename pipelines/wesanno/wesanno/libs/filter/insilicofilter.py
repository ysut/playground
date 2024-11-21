import pandas as pd

from logging import getLogger
logger = getLogger(__name__)


class InSilicoFilter:
    def __init__(self, configs: dict):
        self.configs = configs

    
    def __exclude_low_cadd(self, df: pd.DataFrame) -> pd.DataFrame:
        cutoff = self.configs['Prediction_tools']['cadd']['deleterious']
        df['CADD_phred'].replace('.', 0, inplace=True)
        df = df.astype({'CADD_phred': 'float'})
        df.loc[
            (df['CADD_phred'] >= cutoff) | (df['CADD_phred'] == 0), 
            'CADD_FILTER'] = 'PASS'

        return df

    def __exclude_am_benign(self, df: pd.DataFrame) -> pd.DataFrame:
        """Exclude variants predicted as benign by AlphaMissense
        """
        df.loc[
            (df['AlphaMissense'] == 'ambigous')
              | (df['AlphaMissense'] == 'pathogenic')
              | (df['AlphaMissense'] == '.'),
              'AlphaMissense_FILTER'] = 'PASS'
        
        return df

    def filter(self, df: pd.DataFrame) -> pd.DataFrame:
        logger.info("Filtering in silico predictions")
        

        return df