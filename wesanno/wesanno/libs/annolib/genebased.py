import pandas as pd
import glob


def _genereate_inh_list(remove_inh: str):
    inh_list = ['monoallelic_autosomal', 'monoallelic_X_hem', 
                'monoallelic_Y_hem', 'monoallelic_X_het', 
                'monoallelic_PAR', 'mitochondrial', 
                'biallelic_autosomal', 'biallelic_PAR']
    inh_list.remove(remove_inh)
    
    return inh_list
    

class GeneBasedAnno:
    def __init__(self, args: dict):
        self._root_path = args['resources']
        self.path_to_hgmd = glob.glob(f'{self._root_path}/HGMD*.tsv.gz')
        self.path_to_ddg2p = glob.glob(f'{self._root_path}/DDG2P*.gz')
        self.path_to_eyeg2p = glob.glob(f'{self._root_path}/EyeG2P*.gz')
        self.path_to_sking2p = glob.glob(f'{self._root_path}/SkinG2P*.gz')
        self.path_to_cancer2p = glob.glob(f'{self._root_path}/CancerG2P*.gz')
        self.path_to_cardiacg2p = glob.glob(f'{self._root_path}/CardiacG2P*.gz')
        self.path_to_skeletalg2p = glob.glob(f'{self._root_path}/SkeletalG2P*.gz')

        self.not_ad = _genereate_inh_list('monoallelic_autosomal')
        self.not_ar = _genereate_inh_list('biallelic_autosomal')
        self.not_xl = _genereate_inh_list('monoallelic_X_hem')
        self.not_xl.remove('monoallelic_X_het')
        self.not_yl = _genereate_inh_list('monoallelic_Y_hem')


    def anno_hgmd(self, df: pd.DataFrame) -> pd.DataFrame:
        hgmd = pd.read_table(self.path_to_hgmd[0], header=0, dtype=str)
        hgmd = hgmd[
            ['gene', 'altsymbol', 'refseq', 
             'expected_inheritance', 'hgncID', 'omimid', 'DM']
             ]
        hgmd = hgmd.astype({'DM': 'float64'})
        df = pd.merge(
            df, hgmd, left_on='Gene.refGene', right_on='gene', how='left'
            )
        df = df.drop(columns=['gene'])

        return df
    

    def dm_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        df.loc[df['DM'] != 0, 'non_zero_DM'] = 'PASS'

        return df
    
    
    def anno_dcpr(self, df: pd.DataFrame) -> pd.DataFrame:
        usecols = ['disease name', 'allelic requirement', 'hgnc id']
        panels = ['DD', 'Eye', 'Skin', 'Cancer', 'Cardiac', 'Skeletal']

        dd = pd.read_csv(
            self.path_to_ddg2p[0], header=0, dtype=str, usecols=usecols
            )
        eye = pd.read_csv(
            self.path_to_eyeg2p[0], header=0, dtype=str, usecols=usecols
            )
        ski = pd.read_csv(
            self.path_to_sking2p[0], header=0, dtype=str, usecols=usecols
            )
        can = pd.read_csv(
            self.path_to_cancer2p[0], header=0, dtype=str, usecols=usecols
            )
        car = pd.read_csv(
            self.path_to_cardiacg2p[0], header=0, dtype=str, usecols=usecols
            )
        ske = pd.read_csv(
            self.path_to_skeletalg2p[0], header=0, dtype=str, usecols=usecols
            )

        for df_g2p, panel in zip([dd, eye, ski, can, car, ske], panels):
            df_g2p.rename(columns={
                'gene symbol': 'gene', 'hgnc id': 'hgncID', 
                'disease name': f"{panel}_disease_name",
                'allelic requirement': f"{panel}_allelic_requirment"}, 
                inplace=True
                )
            df_g2p.fillna('.', inplace=True)
            df_g2p = df_g2p.groupby('hgncID').agg(';'.join).reset_index()
            df = pd.merge(df, df_g2p, on='hgncID', how='left')
        
        df['G2P_moi_summary'] = \
            df['DD_allelic_requirment'].fillna('.') + '/' + \
            df['Eye_allelic_requirment'].fillna('.') + '/' + \
            df['Skin_allelic_requirment'].fillna('.') + '/' + \
            df['Cancer_allelic_requirment'].fillna('.') + '/' + \
            df['Cardiac_allelic_requirment'].fillna('.') + '/' + \
            df['Skeletal_allelic_requirment'].fillna('.')
        
        return df


    def anno_gnomad(self, df: pd.DataFrame) -> pd.DataFrame:
        pass


    def summarize_moi(self, row):
        inheritance = row['expected_inheritance']
        alleleinfo = row['G2P_moi_summary']
        
        if inheritance == 'AD':
            return 'AD' if not any(word in alleleinfo for word in self.not_ad) else '.'
        elif inheritance == 'AR':
            return 'AR' if not any(word in alleleinfo for word in self.not_ar) else '.'
        elif inheritance in ['XLD', 'XLR']:
            return 'XL' if not any(word in alleleinfo for word in self.not_xl) else '.'
        elif inheritance == 'YL':
            return 'YL' if not any(word in alleleinfo for word in self.not_yl) else '.'
        else:
            return '.'
    
    

    def match_g2p_phenotypes(self, df):
        pass
