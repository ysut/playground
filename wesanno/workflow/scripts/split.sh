#!/bin/bash

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

INPUT_FILE="$1"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File not found - $INPUT_FILE"
    exit 1
fi

# Split the input file into 4 parts
split -n 4 -d "$INPUT_FILE" "${INPUT_FILE}.part"