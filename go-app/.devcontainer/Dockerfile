# Use official Golang image
FROM golang:1.24.0

# Install some tools
RUN apt-get update && apt-get install -y make git

# 作業ディレクトリ作成 (VS Codeが自動で /workspaces/<プロジェクト名> にマウントしてくれる)
WORKDIR /workspace