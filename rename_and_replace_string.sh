#!/bin/bash
#SBATCH --job-name=mOcc_rename
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=10

# Load command-line arguments
BASE_DIR=$1
OLD_STRING=$2
NEW_STRING=$3

# Check if all arguments are provided
if [[ -z "$BASE_DIR" || -z "$OLD_STRING" || -z "$NEW_STRING" ]]; then
    echo "Usage: $0 <BASE_DIR> <OLD_STRING> <NEW_STRING>"
    exit 1
fi

# Uncompress files in the new directory
gunzip "$BASE_DIR"/*.out.gz
gunzip "$BASE_DIR"/*.align.gz
gunzip "$BASE_DIR"/*.summary.gz
gunzip "$BASE_DIR"/*.divsum.gz

# Function to replace old string within the content of files
replace_in_files() {
    local old_string="$1"
    local new_string="$2"
    find "$BASE_DIR" -type f -exec sed -i "s/$old_string/$new_string/g" {} \;
}

# Replace string in filenames and folder names
find "$BASE_DIR" -depth | while read -r item; do
    # If the old string is found in the name
    if [[ "$(basename "$item")" == *"$OLD_STRING"* ]]; then
        # Replace the old string with the new string in the name
        new_item=$(echo "$item" | sed "s/$OLD_STRING/$NEW_STRING/g")
        
        # Rename the file or directory
        mv "$item" "$new_item"
    fi
done

# Replace old string with new string in the content of all files
replace_in_files "$OLD_STRING" "$NEW_STRING"

# Compress files again
pigz -p 10 "$BASE_DIR"/*.out
pigz -p 10 "$BASE_DIR"/*.align
pigz -p 10 "$BASE_DIR"/*.summary
pigz -p 10 "$BASE_DIR"/*.divsum
