#!/bin/bash
# Script for checking license at begin of files for bash.
# Run with:
#   chmod u+x check_for_license.sh
#   ./check_for_license.sh

# Change directory to the root of the Git repository
cd "$(git rev-parse --show-cdup)"

# Ignore list (directories or subdirectories to ignore)
ignore_list=("cpp/memilio/ad" "_skbuild" "build" "thirdparty", ".*")
print_diff=false

current_year=$(date +'%Y')

# Define the multiline pattern for c++ files with regex for author and contact lines
cpp_pattern=(
    "/\*"
    "\* Copyright \(C\) 2020-$current_year MEmilio"
    "\*"
    "\* Authors:( .+( <[^@<>]+@[^@<>]+>)?(, .+( <[^@<>]+@[^@<>]+>)?)*)?"
    "\*"
    "\* Contact:( .+( <[^@<>]+@[^@<>]+>)?(, .+( <[^@<>]+@[^@<>]+>)?)*)?"
    "\*"
    "\* Licensed under the Apache License, Version 2\.0 \(the \"License\"\);"
    "\* you may not use this file except in compliance with the License\."
    "\* You may obtain a copy of the License at"
    "\*"
    "\*     http://www\.apache\.org/licenses/LICENSE-2\.0"
    "\*"
    "\* Unless required by applicable law or agreed to in writing, software"
    "\* distributed under the License is distributed on an \"AS IS\" BASIS,"
    "\* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied\."
    "\* See the License for the specific language governing permissions and"
    "\* limitations under the License\."
    "\*/"
)

# Define the multiline pattern for python files with regex for author and contact lines
py_pattern=(
    "#+"
    "# Copyright \(C\) 2020-$current_year MEmilio"
    "#"
    "# Authors:( .+( <[^@<>]+@[^@<>]+>)?(, .+( <[^@<>]+@[^@<>]+>)?)*)?"
    "#"
    "# Contact:( .+( <[^@<>]+@[^@<>]+>)?(, .+( <[^@<>]+@[^@<>]+>)?)*)?"
    "#"
    "# Licensed under the Apache License, Version 2\.0 \(the \"License\"\);"
    "# you may not use this file except in compliance with the License\."
    "# You may obtain a copy of the License at"
    "#"
    "#     http://www\.apache\.org/licenses/LICENSE-2\.0"
    "#"
    "# Unless required by applicable law or agreed to in writing, software"
    "# distributed under the License is distributed on an \"AS IS\" BASIS,"
    "# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied\."
    "# See the License for the specific language governing permissions and"
    "# limitations under the License\."
    "#+"
)

# Function to check the first few lines of a file
check_lines() {
    local file=$1
    shift
    local patterns=("$@")
    local line_number=1

    for pattern in "${patterns[@]}"; do
        line=$(sed -n "${line_number}p" "$file")
        if [[ ! "$line" =~ $pattern ]]; then
            echo "$file"
            if $print_diff; then
                diff <(echo "$line") <(echo "$pattern")
            fi
            return
        fi
        ((line_number++))
    done
}

# Find all .cpp, .h, and .py files
files=$(find . \( -name "*.cpp" -o -name "*.h" -o -name "*.py" \))

# Filter out files in the ignore list
filtered_files=()
for file in $files; do
    ignore=false
    for ignore_dir in "${ignore_list[@]}"; do
        if [[ $file == */$ignore_dir/* ]]; then
            ignore=true
            break
        fi
    done
    if ! $ignore; then
        filtered_files+=("$file")
    fi
done

# Check the filtered files
for file in "${filtered_files[@]}"; do
    if [[ $file == *.cpp || $file == *.h ]]; then
        check_lines "$file" "${cpp_pattern[@]}"
    elif [[ $file == *.py ]]; then
        check_lines "$file" "${py_pattern[@]}"
    fi
done