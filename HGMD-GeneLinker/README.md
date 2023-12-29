## NOTE
Linuxだと /home/outputの権限を整えないといけないので，/Docer/app/Dockerfileを修正する必要がある．

## Description
This tiny program is to summrize HGMD data and generate tsv.gz.  
The output file can provide gene-based information such as numbers of DM class variants each gene.

## Usage
1. Download hgmd dump files from HGMD official site (required Professional license)
2. Setup .env file
   Five variables are needed to defined.
   
   1. HGMD version
   2. The path to a directory including HGMD dump files (e.g. hgmd_pro_2023.3.dump.gz)
   3. The path to a directory that binds /var/lib/mysql (For MariaDB data storage)
   4. The path to a directory to save summrized file
   5. The path to a directory to store the session file for phpMyAdmin
   
3. Build images, run containers, and login container by folowing codes
```
$ docker compose up -d
$ docker compose exec -it app bash
```
In this step, log file are created in your host directory defined by ${BIND_DIR}.

4. Activate conda environment (named maria) and run tasks
```
(base) # conda activate maria
(maria) # task setup
(maria) # task summarize
```


#### phpMyAdmin
The phpMyAdmin container will run. You can access `http://localhost:3000/` from host using web browser.

