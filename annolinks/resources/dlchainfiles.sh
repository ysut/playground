#!bin/bash

set -e

readonly BASE_URI="ftp://hgdownload.soe.ucsc.edu/goldenPath"
readonly hg19tohg38="hg19ToHg38.over.chain.gz"
readonly hg38tohg19="hg38ToHg19.over.chain.gz"

# Download function using curl
function download_with_curl() {
  echo "Downloading with curl: $1"
  filename="$1"
  assembly="${filename:0:4}"
  curl -R -o "$1" "${BASE_URI}"/"${assembly}"/liftOver/"$1"
}

# Download function using wget
function download_with_wget() {
  echo "Downloading with wget: $1"
  filename="$1"
  assembly="${filename:0:4}"
  wget --timestamping "${BASE_URI}"/"${assembly}"/liftOver/"$1" -O "$1"
}

function download_file() {
  if command -v curl > /dev/null; then
    download_with_curl "$1"
  elif command -v wget > /dev/null; then
    download_with_wget "$1"
  else
    echo "Neither curl nor wget is installed."
    exit 1
  fi
}

function download_chainfiles() {
  download_file "${hg19tohg38}"
  download_file "${hg38tohg19}"
}

function main() {
  download_chainfiles
}

main
