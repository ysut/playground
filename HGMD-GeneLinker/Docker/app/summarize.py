import sys

import pandas as pd
import pymysql


## Set HGMD Version for output
hgmd_ver: str = sys.argv[1]
output_path = f"{sys.argv[2]}/HGMD_GeneBasedInfo_{hgmd_ver}.tsv.gz"


## Connection
database = 'hgmd_pro'
mariadb_params = {
    'host': 'db', 'user': 'root', 'password': 'pwd',
    'port': 3306, 'database': database, 
    'cursorclass': pymysql.cursors.DictCursor
    }
connection = pymysql.connect(**mariadb_params)


## Create tags DataFrame
sql = f"SELECT gene, tag FROM allmut;"
with connection.cursor() as cursor:
    cursor.execute(sql)
    fetched_allmut = cursor.fetchall()

df_allmut = pd.DataFrame(fetched_allmut)

df_tags = df_allmut.pivot_table(
    index='gene', columns='tag', aggfunc='size', fill_value=0
    )
df_tags = df_tags.reset_index()
df_tags.columns.name = None


## Create all genes DataFrame
sql = f"SELECT gene, altsymbol, refseq, expected_inheritance, hgncID, omimid FROM allgenes;"

with connection.cursor() as cursor:
    cursor.execute(sql)
    fetched_allgenes = cursor.fetchall()

df_genes = pd.DataFrame(fetched_allgenes)

## Left join
df = pd.merge(left=df_genes, right=df_tags, how='left', on='gene')


## Save as pickle
df.to_csv(
    output_path, index=False, sep='\t', encoding='utf-8', compression='gzip'
    )