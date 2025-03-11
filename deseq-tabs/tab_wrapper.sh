#!/bin/bash
# Change to your target folder
#cd data

# List all .tab files, wrap each in quotes, and join with commas
files=$(ls *.tab | jq -R . | jq -s .)

# Save the JSON array to tabFiles.json in the data folder
echo "$files" > tabFiles.json

# Optionally, print the result
cat tabFiles.json
