# Description
To annotate hyperlinks such as HGMD, UCSC genome browser, DECIPHER for each variants.

# Required input file columns
CHROM, POS, REF, ALT and Gene.  
POS column is hg38.

# Usage

Simple example:
```Shell
$ annolinksver=1.2
$ cd /path/to/your/directory  # Move to the directory where the input files are located.
$ ls
Sample.xlsx

$ docker run --rm -v $(pwd):/input utsuno/annolinks:${annolinksver} \
  python -m annolinks \
  --input /input/Sample.xlsx \
  --gene-col Gene.refGene
```

## Arguments and options
| Arguments | Required | Description | Default |
| ---  | --- | ---- | --- |
| `--input` / `-I`   | True | Input file path   | None         |
| `--gene-col` / `-G` | True | Gene Symbols column name | None|
| `--output` / `-O`  | False | Output file path  | /SAME_AS_INPUT/INPUT_FILE_hyperlinked.xlsx|
| `--windows` | False | Using output Excel on windows | False|
| `--franklin-page` | False | Choice one from <br>`assessment-tools`,<br>`variant-interpretation`,<br>`publications`,<br>`gene-assessment`,<br>`conditions`,<br>`clinical-evidence`,<br>`community`,<br>and `classification-demo-app`| `assessment-tools` |
| `--spliceai-raw` | False | Switch between mask score and raw score | False (masked score)|
| `--spliceai-dist` | False | SpliceAI max distance | `10000` |
| `--skip-sheets` | False | Specify sheets not to be annoteted.<br>| None|
| `--skip-sites` | False |||
| ...|...|...|...|


# Usage in Japanese
## NOTE
- Excelファイルに対応しています．
- 少なくともCHROM，POS，REF，ALT，Gene（列名は任意）の列が必要です．
- POS列は，hg38であると特別なオプションは必要なくすべての機能を実行できます．
- SpliceAI lookupのURLをアノテーションする際，アウトプットされたExcelファイルをWindowsで扱う場合は，`--windows`オプションを付けてください．



