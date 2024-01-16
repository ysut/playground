## Resources for wesanno
1. Gene-based HGMD information file
2. Gene2Phenotype CSVs

      

## How do I prepare these resources?
1. Gene-based HGMD information file  
    1.1 Using my tiny program.  
    https://github.com/ysut/playground/tree/main/HGMD-GeneLinker

    1.2 Copy a tsv.gz to resources directory

2. Gene2Phenotype CSVs  
    2.1 Install `wget` or `curl`.  
    2.2 Using a shell script (dlg2p.sh).  
    ```Shell
    $ cd resources
    $ ./dlg2p.sh
    ```
