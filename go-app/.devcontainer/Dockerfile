# ベースイメージに公式Goイメージを使う
FROM golang:1.23.4

# 必要に応じてツール追加
# RUN apt-get update && apt-get install -y make git ...

# 作業ディレクトリ作成 (VS Codeが自動で /workspaces/<プロジェクト名> にマウントしてくれる)
WORKDIR /workspace