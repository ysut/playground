# Description
To annotate hyperlinks such as HGMD, UCSC genome browser, DECIPHER for each variants.

# Required input file columns
CHROM, POS, REF, ALT and Gene.  
POS column is hg38.

# Usage

Simple example on Windows PC:
```Shell
$ cd /path/to/your/directory  # Move to the directory where the input files are located.
$ ls 
sample.xlsx

$ docker run --rm -v $(pwd):/input utsuno/annolinks:latest \
  python -m annolinks \
  --input /input/sample.xlsx \    # Please change "sample.xlsx" to your excel file name. The "/input/" path isn't needed to change. 
  --gene-col Gene.refGene \
  --windows
```

## Arguments and options
| Arguments | Required | Description | Default |
| ---  | --- | ---- | --- |
| `--input` / `-I`   | True | Input file path   | None         |
| `--gene-col` / `-G` | True | Gene Symbols column name | None |
| `--output` / `-O`  | False | Output file path  | /SAME_AS_INPUT/INPUT_FILE_hyperlinked.xlsx|
| `--windows` | False | Using output Excel on windows | None |
| `--franklin-page` | False | Choice one from <br>`assessment-tools`,<br>`variant-interpretation`,<br>`publications`,<br>`gene-assessment`,<br>`conditions`,<br>`clinical-evidence`,<br>`community`,<br>and `classification-demo-app`| `assessment-tools` |
| `--spliceai-raw` | False | Switch between mask score and raw score | False (masked score)|
| `--spliceai-dist` | False | SpliceAI max distance | `10000` |
| `--skip-sheets` | False | Specify sheets not to be annoteted.| None (All sheets will be annotated.)|
| `--skip-sites` | False | Specify sites not to be inserted hyperlinks.|None (Hyperlinks to all sites will be inserted.)|
| ...|...|...|...|


# Usage in Japanese
## 必要最低限の情報
- .xlsxファイルに対応しています．
- Excelファイルには，少なくとも次の5列が必要です．
  "CHROM"，"POS"，"REF"，"ALT"，と「遺伝子名が書いてある列（列名は任意）」が必要です．そして
  「遺伝子名が書いてある列」は，`--gene-col1`で必ず指定する必要があります．
  例えば，`--gene-col Gene.refGene`など．
- ハイパーリンクを挿入したExcelファイルをWindowsで扱う場合は，`--windows`オプションを付けてください．
  MacOSやLinuxでExcelファイルを見る場合は必要ありません．
- 遺伝子名の列に「MARCH1」とか「SEPT7」とか日付と間違えるものが含まれているとエラーが出ます．
