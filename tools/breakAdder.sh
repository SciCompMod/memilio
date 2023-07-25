cd $(git rev-parse --show-cdup) # Change to the root of the repository
git ls-files -z | # List all files in the repository
while IFS= read -rd '' f; do # Read each file in the repository
     if [[ $f != *.(png|md|h5|txt) ]]; then
        tail -c1 < "$f" | read -r _ || echo >> "$f"; # Add a newline if the file has no end of line break in the last line
    fi
done
