#! /bin/bash
set -e

for i in {5..7}; do
    echo "Running Analysis ${i}"
    python "/home/mcaskey/10XvParse/Scripts/analysis${i}.py"
    echo "Analysis ${i} completed."
done

echo "All scripts completed."