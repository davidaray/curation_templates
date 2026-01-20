#!/bin/bash
#SBATCH --job-name=mVel_rename
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1



# Set the base directory, old string, and new string
BASE_DIR="/lustre/scratch/daray/bat1k_TE_analyses/hite/mMyoVel_1"
OLD_STRING="mMyoVel1"
NEW_STRING="mMyoVel1.pri"

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

