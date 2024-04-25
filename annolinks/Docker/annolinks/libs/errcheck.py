import dataclasses
import pandas as pd
from logging import getLogger

logger = getLogger(__name__)

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

    def has_specified_columns(self) -> bool:
        specified_cols = [
            'HGMD', 'DECIPHER', 'SpliceAI_Lookup', 'UCSC', 'Franklin'
            ]
        
        result = {}
        for sheet in self.anno_sheets:
            must_change_cols = []
            for col in specified_cols:
                if col in self.dfs.sheets[sheet].columns:
                    must_change_cols.append(col)
            result[sheet] = must_change_cols

        if sum([len(v) for v in result.values()]):
            logger.error(f'The following column names must be changed: '
                         f'{result}')
            return True
        else:
            logger.info('OK! All specified columns are not found.')
            return False

    def has_required_columns(self) -> bool:
        required_cols = ['CHROM', 'POS', 'REF', self.args['alt_col']]
        result = {}        
        for sheet in self.anno_sheets:
            must_cols = []
            for col in required_cols:
                if col not in self.dfs.sheets[sheet].columns:
                    must_cols.append(col)
            result[sheet] = must_cols

        if sum([len(v) for v in result.values()]):
            logger.error(f'The following column must be required to annotate: '
                         f'{result}')
            return False
        else:
            logger.info('OK! All required columns are found.')
            return True
            
    # def __has_liftover_columns(self) -> bool:
    #     if ((self.args['assembly'] == 'hg19') 
    #         | (self.args['assembly'] == 'GRCh37')):
    #         if not 'DECIPHER' in self.args['--skip-sites']:
    
    def __count_links(self) -> int:
        anno_sites_num = 5 - len(self.skip_sites)        
        total = 0
        for sheet in self.anno_sheets:
            total = total + self.dfs.sheets[sheet].shape[0] * anno_sites_num
        return total

    def is_within_limit(self) -> bool:
        """
        Check if the number of total links in the excel workbook exceeds the limit.
        """
        total = self.__count_links()
        if total <= 65530:
            logger.info(f"The total number of hyperlinks is {total} (≤ 65530).")
            return True
        else:
            logger.error("The total number of hyperlinks exceeds the limit (65530): "
                         f"The total number of hyperlinks is {total}. "
                         "Please reduce the input variants or"
                         "consider using the '--skip-sheets' or "
                         "'--skip-sites' options.")
            return False

    # Check if the cell value of input excel file is datetime data type
    def exists_datetime(self) -> bool:
        for sheet in self.anno_sheets:
            if self.dfs.sheets[sheet].dtypes.any() == 'datetime64[ns]':
                logger.error(f"Sheet '{sheet}' contains datetime data type.")
                return False
        logger.info("OK! All gene symbol columns do not contain datetime data type.")
        return True
