#!/bin/bash

echo "STEP 2"

input="$1"

awk -F '\t' '{print $1":"$2"-"$4"-"$5}' "$1"
