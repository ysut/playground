#!/bin/bash

# This script is used to setup the HGMD data in the MySQL database.
# After login the docker ocntainer, run this script to setup the HGMD data.

set -eu

SCRIPT_DIR=$(cd $(dirname $0); pwd)

# cd
cd "${HGMD_DUMPS_DIR}"

# Decompress the HGMD data files
echo "Decompressing HGMD data files..."
if [ -f "hgmd_snp-${HGMD_VERSION}.dump" ]; then
  echo "[INFO] hgmd_snp-${HGMD_VERSION}.dump has already been decompressed."
else
  gunzip hgmd_snp-${HGMD_VERSION}.dump.gz
fi
if [ -f "hgmd_pro-${HGMD_VERSION}.dump" ]; then
  echo "[INFO] hgmd_pro-${HGMD_VERSION}.dump has already been decompressed."
else
  gunzip hgmd_pro-${HGMD_VERSION}.dump.gz
fi
if [ -f "hgmd_phenbase-${HGMD_VERSION}.dump" ]; then
  echo "[INFO] hgmd_phenbase-${HGMD_VERSION}.dump has already been decompressed."
else
  gunzip hgmd_phenbase-${HGMD_VERSION}.dump.gz
fi
if [ -f "hgmd_views-${HGMD_VERSION}.dump" ]; then
  echo "[INFO] hgmd_views-${HGMD_VERSION}.dump has already been decompressed."
else
  gunzip hgmd_views-${HGMD_VERSION}.dump.gz
fi

# Reset database for the HGMD data
mysql --defaults-extra-file="${SQL_CNF}" \
  --execute "DROP DATABASE IF EXISTS hgmd_pro;"
mysql --defaults-extra-file="${SQL_CNF}" \
  --execute "DROP DATABASE IF EXISTS hgmd_snp;"
mysql --defaults-extra-file="${SQL_CNF}" \
  --execute "DROP DATABASE IF EXISTS hgmd_phenbase;"
mysql --defaults-extra-file="${SQL_CNF}" \
  --execute "DROP DATABASE IF EXISTS hgmd_views;"

# Create database for the HGMD data
mysqladmin --defaults-extra-file="${SQL_CNF}" create hgmd_snp
mysqladmin --defaults-extra-file="${SQL_CNF}" create hgmd_pro
mysqladmin --defaults-extra-file="${SQL_CNF}" create hgmd_phenbase
mysqladmin --defaults-extra-file="${SQL_CNF}" create hgmd_views

# Restore the HGMD data
echo "Restoring HGMD data..."

echo "Restoring hgmd_snp..."
mysql --defaults-extra-file="${SQL_CNF}" hgmd_snp \
  < hgmd_snp-${HGMD_VERSION}.dump

echo "Restoring hgmd_pro..."
mysql --defaults-extra-file="${SQL_CNF}" hgmd_pro \
  < hgmd_pro-${HGMD_VERSION}.dump

echo "Restoring hgmd_phenbase..."
mysql --defaults-extra-file="${SQL_CNF}" hgmd_phenbase \
  < hgmd_phenbase-${HGMD_VERSION}.dump

echo "Restoring hgmd_views..."
mysql --defaults-extra-file="${SQL_CNF}" hgmd_views \
  < hgmd_views-${HGMD_VERSION}.dump

cd $SCRIPT_DIR
