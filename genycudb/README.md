
# Description

For YCU lab.
To summarize sample information as a SQLite database.


# Usage (in Japanese)
### 1. Copy three Excel files  
rawdataディレクトリに以下の3つのエクセルファイルをコピーする．

  1. 古いREP用のエクセルファイル
  2. 現在更新中のREP用のエクセルファイル
  3. 郵送済みリストのエクセルファイル

### 2. 
Setup the environment using the following commands in a conda-enabled setting.
```Shell
$ conda create -n ycudb python=3.9 -y && conda activate ycudb
$ conda install -y pandas 
```

### 3. RUN
```Shell
(ycudb) $ python main.py
```
dbディレクトリに SQLite database ができている．

#### NOTE
エクセルファイル内の入力をうまく整理できなくて途中で止まるかもしれない．