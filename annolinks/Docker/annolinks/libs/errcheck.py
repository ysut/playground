import dataclasses
import pandas as pd
from logging import getLogger

logger = getLogger(__name__)

# Check required columns in input excel file
# Minimum required columns are 'CHROM', 'POS', 'REF', 'ALT'

class ErrorCheck:
    def __init__(
            self,
            dfs: dataclasses.dataclass,
            anno_sheets: list,
            skip_sites: list,
            args: dict
            ):
        self.dfs = dfs
        self.anno_sheets = anno_sheets
        self.skip_sites = skip_sites
        self.args = args

    def __has_required_columns(self, df: pd.DataFrame) -> bool:
        required_cols = ['CHROM', 'POS', 'REF', 'ALT']
        for col in required_cols:
            if col not in df.columns:
                return False
        return True
    
    def __has_gene_column(self, df: pd.DataFrame) -> bool:
        if self.args['gene_col'] in df.columns:
            return True
        else:
            return False
    
    # If gene column is not in the input excel file,
    # HGMD, 
    def __adjust_skip_sites(self) -> list:
        pass


    # def __has_liftover_columns(self) -> bool:
    #     if ((self.args['assembly'] == 'hg19') 
    #         | (self.args['assembly'] == 'GRCh37')):
    #         if not 'DECIPHER' in self.args['--skip-sites']:
    
    def __count_links(self) -> int:
        """
        Check if the number of total links in the excel workbook exceeds the limit.
        """
        anno_sites_num = 5 - len(self.skip_sites)
        
        total = 0
        for sheet in self.anno_sheets:
            total = total + self.dfs.sheets[sheet].shape[0] * anno_sites_num

        return total

    def is_within_limit(self) -> bool:
        total = self.__count_links()
        if total <= 65530:
            logger.info(f"Total number of hyperlinks are {total}.")
            return True
        else:
            logger.error("The number of total links exceeds the limit (65530): "
                         f"Total number of hyperlinks are {total}."
                         "Please reduce the number of input variants or"
                         "consider using the '--skip-sheets' or "
                         "'--skip-sites' options.")
            return False
