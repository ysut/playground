# Description
To annotate hyperlinks such as HGMD, UCSC genome browser, DECIPHER, and SpliceAI lookup for each variants.

# Required input file columns
"CHROM", "POS", "REF", "ALT" and gene symbol column. 

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
| `--input` / `-I`   | Required | Input file path   | None         |
| `--gene-col` / `-G` | Required | Gene Symbols column name | None |
| `--pos19` | No | The column name of variant position on hg19 or GRCh37 | `POS`|
| `--pos38` | No | The column name of variant position on hg38 or GRCh38 | None|
| `--assembly` / `-A`   | No | Select assembly version | `hg19`        |
| `--output` / `-O`  | No | Output file path  | /SAME_AS_INPUT/INPUT_FILE_hyperlinked.xlsx|
| `--windows` | No | Using output Excel on windows | None |
| `--franklin-page` | No | Choice one from <br>`assessment-tools`,<br>`variant-interpretation`,<br>`publications`,<br>`gene-assessment`,<br>`conditions`,<br>`clinical-evidence`,<br>`community`,<br>and `classification-demo-app`| `assessment-tools` |
| `--spliceai-raw` | No | Switch between mask score and raw score | False (masked score)|
| `--spliceai-dist` | No | SpliceAI max distance | `10000` |
| `--skip-sheets` | No | Specify sheets not to be annoteted.| None (All sheets will be annotated.)|
| `--skip-sites` | No | Specify sites not to be inserted hyperlinks.|None (Hyperlinks to all sites will be inserted.)|
| `--ucsc-width` | No | Number of bases on both sides of the variant in the UCSC Genome Browser | `45` |
| `--split-alt` | No | If ALT has multiple alternatives, separate each variant into its own row.<br> Please note that this will increase the number of rows. | False |
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
