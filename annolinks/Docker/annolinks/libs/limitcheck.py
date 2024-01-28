import dataclasses
from logging import getLogger

logger = getLogger(__name__)

class LimitCheck:
    def __init__(
            self,
            dfs: dataclasses.dataclass,
            anno_sheets: list,
            skip_sites: list
            ):
        self.dfs = dfs
        self.anno_sheets = anno_sheets
        self.skip_sites = skip_sites

    def count_links(self) -> int:
        """
        Check if the number of total links in the excel workbook exceeds the limit.
        """
        anno_sites_num = 5 - len(self.skip_sites)
        
        total = 0
        for sheet in self.anno_sheets:
            total = total + self.dfs.sheets[sheet].shape[0] * anno_sites_num

        return total

    def is_within_limit(self) -> bool:
        total = self.count_links()
        if total <= 65530:
            return True
        else:
            logger.error("The number of total links exceeds the limit (65530): "
                         f"Total number of hyperlinks are {total}."
                         "Please reduce the number of input variants or"
                         "consider using the '--skip-sheets' or "
                         "'--skip-sites' options.")
            return False
