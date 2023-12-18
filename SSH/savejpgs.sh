#!/bin/bash

set -eu

readonly OUTPUT_DIR=$(realpath ${1})
readonly BASE_URL="https://s3-ap-northeast-1.amazonaws.com/cdn.jrha.or.jp/img"

# Make a directory to save the jpgs
mkdir -p "${OUTPUT_DIR}" && cd "$_"

# Set the url array
url_array="${BASE_URL}/2021/[001-248].jpg ${BASE_URL}/2021/[301-541].jpg \
           ${BASE_URL}/2020/[001-258].jpg ${BASE_URL}/2020/[301-539].jpg \
           ${BASE_URL}/2019/[001-246].jpg ${BASE_URL}/2019/[301-526].jpg \
           ${BASE_URL}/2018/[001-243].jpg ${BASE_URL}/2018/[301-539].jpg \
           ${BASE_URL}/2017/[001-247].jpg ${BASE_URL}/2017/[301-528].jpg \
           ${BASE_URL}/2016/[001-248].jpg ${BASE_URL}/2016/[301-543].jpg \
           ${BASE_URL}/2015/[001-244].jpg ${BASE_URL}/2015/[301-541].jpg \
           ${BASE_URL}/2014/[001-260].jpg ${BASE_URL}/2014/[301-531].jpg \
           ${BASE_URL}/2013/[001-261].jpg ${BASE_URL}/2013/[301-529].jpg \
           ${BASE_URL}/2012/[001-250].jpg ${BASE_URL}/2012/[301-527].jpg \
           ${BASE_URL}/2011/[001-240].jpg ${BASE_URL}/2011/[301-528].jpg \
           ${BASE_URL}/2010/[001-220].jpg ${BASE_URL}/2010/[301-515].jpg \
           ${BASE_URL}/2009/[001-161].jpg ${BASE_URL}/2009/[201-536].jpg \
           ${BASE_URL}/2008/[001-163].jpg ${BASE_URL}/2008/[201-536].jpg"

# Download the jpgs in the directory of each year
for url in ${url_array}; do
  mkdir -p "${OUTPUT_DIR}/$(basename $(dirname ${url}))" && cd "$_"
  curl \
    --silent \
    --show-error \
    --location \
    --fail \
    --retry 3 \
    -O "${url}"
  cd "${OUTPUT_DIR}"
  sleep $(($RANDOM % 20))
done

mkdir -p "${OUTPUT_DIR}"/download.done
