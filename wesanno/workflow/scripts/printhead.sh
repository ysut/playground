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

# Print the first 10 lines of the input file
head -n 10 "$INPUT_FILE"
