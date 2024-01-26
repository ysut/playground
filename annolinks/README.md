# Description
To annotate hyperlinks such as HGMD, UCSC genome browser, DECIPHER for each variants.

# Required input file columns
CHROM, POS, REF, ALT and Gene.
POS column is hg38.

# Usage

Example:
```
$ cd /path/to/your/directory
$ ls
Sample.xlsx

$ docker run --rm -v $(pwd):/input utsuno/repo:annolinks \
  python -m annolinks \
  --input /input/Sample.xlsx \
  --gene-col Gene.refGene
```
