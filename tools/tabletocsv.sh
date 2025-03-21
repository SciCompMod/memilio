#!/bin/bash

# Use this script to convert the table output of MEmilio into a csv-file.
# Usage: Save the table output as a file, then call this script with
#
#    bash table_file csv_file
#
# where you also need to add the correct file endings.

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 input_file output_file"
  exit 1
fi

sed '/^$/d' "$1" |        # Remove empty lines
sed 's/^[ \t]*//' |                # Remove leading spaces/tabs
sed 's/[ \t]*$//' |                # Remove trailing spaces/tabs
tr -s '[:blank:]' ';' > "$2"  # Replace multiple spaces/tabs with a single semicolon
