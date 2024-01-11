## Description
This tiny program is to summarize HGMD data and generate tsv.gz.
The output file can provide gene-based information such as the number of DM class variants in each gene.

## Usage
1. Download hgmd dump files from HGMD official site (required Professional license)
2. Move to the `HGMD-GeneLinker` directory.  
   Setup `.env` file. Five variables are needed.
   
   1. HGMD version
   2. The path to a directory including HGMD dump files (e.g. hgmd_pro_2023.3.dump.gz)
   3. The path to a directory that binds /var/lib/mysql (For MariaDB data storage)
   4. The path to a directory to save summrized file
   5. The path to a directory to store the session file for phpMyAdmin
   
4. Build images, run containers, and login container by folowing codes.
```
$ docker compose up -d
$ docker compose exec app bash
```

4. Activate conda environment (named maria) and run tasks  
In this step, a log file (e.g. log-HGMD.2023.3.txt) are created in the host directory defined by ${BIND_DIR}.
```
(base) # task
```
It will take a long time⏳.  
Finaly, a gene-based summarized tsv (e.g. HGMD_GeneBasedInfo_2023.3.tsv.gz) will be created.


## Output format
Columns are:
gene, altsymbol, refseq, expected_inheritance, hgncID, omimid, DFP, DM, DM?, DP, FP, R


## phpMyAdmin
The phpMyAdmin container will run. You can access `http://localhost:3000/` from host using web browser.  
If you can not access phpMyAdmin, you should check the permission of sessions directory.  
Session directory permissions must be 777.

