import os
import collections
from pathlib2 import Path

import pysam
import numpy as np
import pandas as pd

from pandarallel import pandarallel
os.environ['JOBLIB_TEMP_FOLDER'] = '/tmp' 
pandarallel.initialize(progress_bar=True, verbose=1)

class Anno:
    def __init__(self, df: pd.DataFrame, args: dict):
        self.df: pd.DataFrame = df
        self.args: dict = args
        self.resources_dir: Path = args['resources']

    #--- ENST ID ---
    def __anno_enst(self, row) -> str:
        query_region: str = f"chr{row['CHROM']}:{row['Start']}"
        
        pass

    #--- HGMD gene-based info ---
    def __load_hgmd_info(self) -> pd.DataFrame:
        hgmd_resource: str = [
            str(x) for x in self.resources_dir.glob(
                'HGMD/HGMD_gene_based.tsv.gz'
                )
            ][0]
        return pd.read_csv(hgmd_resource, sep='\t')
    


    #--- gnomAD Constraint ---
    def __gnomad_constraint(self, row) -> str:
        pass



    #--- AlphaMissense ---
    def __load_alphamissense(self) -> pysam.TabixFile:
        am_resource: str = [
            str(x) for x in self.resources_dir.glob(
                'AlphaMissense/AlphaMissense_hg19.tsv.gz'
                )
            ][0]
        return pysam.TabixFile(am_resource)
    
    def __fetch_am_info(self, row) -> str:
        tbx_am = self.__load_alphamissense()
        chrom, pos = f"chr{row['CHROM']}", int(row['POS'])
        ref, alt = row['REF'], row['ALT']
        if ((len(ref) == 1) & (len(alt) == 1)):            
            am_rows = tbx_am.fetch(
                chrom, pos-1, pos, parser=pysam.asTuple())
            for am_row in am_rows:
                if am_row[3] == alt:
                    return f'{am_row[8]}:{am_row[9]}'
                else:
                    pass
        else:
            pass # If indel, return missing value
        
    def __separate_am_info(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[df['AM'].isna(), 'AM'] = '.:.'
        df_am = df['AM'].str.split(':', expand=True)
        df_am.columns = ['AM_score', 'AM_phred']
        df = pd.concat([df, df_am], axis=1)
        df = df.drop('AM', axis=1)

        return df

    def __anno_alphamissense(self, df: pd.DataFrame) -> pd.DataFrame:
        df['AM'] = df.progress_apply(self.__fetch_am_info, axis=1)
        df = self.__separate_am_info(df)
 
        return df


    #--- REVEL ---
    def __load_revel(self):
        revel_resource: str = [
            str(x) for x in self.resources_dir.glob(
                'REVEL/REVEL_v1.3_all_chromosomes.tsv.gz'
                )
            ][0]
        self.tbx_revel = pysam.TabixFile(revel_resource)

    def __fetch_revel_scores(
            self, chrom: str, pos: int, ref: str, alt: str, enst: str) -> str:
        revel_rows = self.tbx_revel.fetch(
            chrom, pos-1, pos, parser=pysam.asTuple()
            )
        
        for revel_row in revel_rows:
            if revel_row[4] == alt:
                if revel_row[7] == enst:
                    # return revel_row[4]
                    pass
                else:
                    pass
            else:
                pass

    def __anno_revel(self, row):
        score = self.__fetch_revel_scores(
            chroom=row['CHROM'], pos=int(row['POS']), 
            ref=row['REF'], alt=row['ALT'], enst=row['ENST']
            )



    #--- Trap ---
    def __anno_trap(self):
        pass

    #--- Synvep ---
    def __anno_synvep(self):
        pass
    def __anno_eve(self):
        pass
    def __anno_shine(self):
        pass
    def __anno_maverick(self):
        pass

    #--- Main function ---
    def anno_scores(self) -> pd.DataFrame:
        if self.args['no_alphamissense']:
            pass
        else:
            self.df = self.__anno_alphamissense(self.df)

        # if self.args['no_revel']:
        #     pass
        # else:
        #     self.__load_revel()
        #     self.df['REVEL'] = self.df.progress_apply(
        #         self.__anno_revel, axis=1)

        # if self.args['no_trap']:
        #     pass
        # else:
        #     pass
        
        return self.df


#==============================================================================#
#==============================================================================#
#==============================================================================#


class MaverickAnnotator:
    def __init__(self, pre_computed_scores: str):
        self.pre_computed_scores = pre_computed_scores
        self.tbx_mav = pysam.TabixFile(self.pre_computed_scores)

    def _fetch_mav_scores(
            self, chrom: str, pos: int, alt: str) -> collections.namedtuple:
        
        if (chrom == 'X') | (chrom == 'Y') | (chrom == 'MT'):
            return '.'
        else:
            mav_rows = self.tbx_mav.fetch(
                chrom, pos-1, pos, parser=pysam.asTuple())
            
            for mav_row in mav_rows:
                if mav_row[4] == alt:
                    MavScores = collections.namedtuple(
                        'MavScores', ['benign', 'dominant', 'recessive'])
                    return MavScores(mav_row[9], mav_row[10], mav_row[11])
                else:
                    pass
            return '.'
            
    def anno_benign_to_df(self, row):
        if ((len(row['REF']) == 1) & (len(row['ALT']) == 1)):
            mav_scores = self._fetch_mav_scores(
                chrom=row['CHROM'], pos=int(row['POS']), alt=row['ALT'])
            # print(f"{row['CHROM']}:{row['POS']}-{row['REF']}-{row['ALT']}")
            # print(mav_scores)
            if mav_scores == '.':
                return mav_scores
            else:
                return mav_scores.benign     
    
    def anno_dominant_to_df(self, row):
        if ((len(row['REF']) == 1) & (len(row['ALT']) == 1)):
            mav_scores = self._fetch_mav_scores(
                chrom=row['CHROM'], pos=int(row['POS']), alt=row['ALT'])
            if mav_scores == '.':
                return mav_scores
            else:
                return mav_scores.dominant

    def anno_recessive_to_df(self, row):
        if ((len(row['REF']) == 1) & (len(row['ALT']) == 1)):
            mav_scores = self._fetch_mav_scores(
                chrom=row['CHROM'], pos=int(row['POS']), alt=row['ALT'])
            if mav_scores == '.':
                return mav_scores
            else:
                return mav_scores.recessive



class TrapAnnotator:
    # chr19   54484   54485   A       T       ENSG00000248385 0.017
    def __init__(self, pre_computed_scores_dir: str) -> None:
        self.pre_computed_scores_dir = pre_computed_scores_dir
        # self.tbx_trap = pysam.TabixFile(self.pre_computed_scores, )

    def liftover_to_hg38(self, row) -> int:
        converter = get_lifter('hg19', 'hg38')
        result = converter.query(row['CHROM'], row['POS'])
        if result:
            return int(result[0][1])
        else:
            return None
        
    def _fetch_trap_scores(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        trap_scores = (f'{self.pre_computed_scores_dir}/'
                       '{chrom}.hg38.tsv.bed.sorted.gz')
        tbx_trap = pysam.TabixFile(trap_scores)

        for trap_row in tbx_trap.fetch(
            f"chr{chrom}", pos-1, pos, parser=pysam.asBed()):
            if trap_row[4] == alt:
                print(trap_row)
                # return trap_row[6]


    def anno_to_df(self, row, pos_hg38_col: str) -> str:
        trap_score = self._fetch_trap_scores(
            row['CHROM'], row['POS'], row['REF'], row['ALT'])
    
        pass


class SynvepAnnotator:
    def __init__(self, db: str) -> None:
        self.db = db
        self.conn = sqlite3.connect(self.db)
        self.cur = self.conn.cursor()

    def _fetch_synvep_scores(
            self, chrom: str, pos: int, ref: str, alt: str, enst:str):
        sql = f"SELECT * \
                FROM VARIANT_SCORE \
                WHERE Ensembl_transcript_ID='{enst}' AND \
                chr='{chrom}' AND \
                pos='{pos}' AND \
                ref='{ref}' AND \
                alt='{alt}';"
        self.cur.execute(sql)
        fetched_data = self.cur.fetchone()
        return fetched_data

    def anno_to_df(self, row, enst_col: str):
        enst = row[enst_col]
        chrom, pos = str(row['CHROM']), str(row['POS'])
        ref, alt = row['REF'], row['REF']
        synvep_scores = self._fetch_synvep_scores(chrom, pos, ref, alt, enst)

        return synvep_scores[8]


### ## ### ## ### ## ### ## ### ## ### ## ### ## ### ## ### ## ### ## ### ## ###
class EveAnnotator:
    def __init__(self) -> None:
        pass

class ShineAnnotaor:
    def __init__(self) -> None:
        pass

### ## ### ## ### ## ### ## ### ## ### ## ### ## ### ## ### ## ### ## ### ## ###
