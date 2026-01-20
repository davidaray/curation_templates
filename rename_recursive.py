import os
import argparse

def RENAME_WITHIN_FILE(FILE_PATH, OLD_NAME, NEW_NAME):
    """Replace occurrences of OLD_NAME with NEW_NAME within the file."""
    try:
        with open(FILE_PATH, 'r') as FILE:
            CONTENT = FILE.read()
        
        NEW_CONTENT = CONTENT.replace(OLD_NAME, NEW_NAME)
        
        if NEW_CONTENT != CONTENT:
            with open(FILE_PATH, 'w') as FILE:
                FILE.write(NEW_CONTENT)
            print(f"Replaced content in: {FILE_PATH}")
    except (UnicodeDecodeError, IOError):
        print(f"Skipping binary or inaccessible file: {FILE_PATH}")

def RENAME_FILES_FOLDERS(ROOT_DIR, OLD_NAME, NEW_NAME):
    """Recursively rename only files, folders, and content containing the old_name."""
    for DIRPATH, DIRNAMES, FILENAMES in os.walk(ROOT_DIR, topdown=False):
        # Filter and rename files that contain the old_name
        for FILENAME in [f for f in FILENAMES if OLD_NAME in f]:
            OLD_FILE_PATH = os.path.join(DIRPATH, FILENAME)
            NEW_FILENAME = FILENAME.replace(OLD_NAME, NEW_NAME)
            NEW_FILE_PATH = os.path.join(DIRPATH, NEW_FILENAME)
            os.rename(OLD_FILE_PATH, NEW_FILE_PATH)
            print(f"Renamed file: {OLD_FILE_PATH} -> {NEW_FILE_PATH}")
            RENAME_WITHIN_FILE(NEW_FILE_PATH, OLD_NAME, NEW_NAME)
        
        # Rename directories that contain the old_name
        for DIRNAME in [d for d in DIRNAMES if OLD_NAME in d]:
            OLD_DIR_PATH = os.path.join(DIRPATH, DIRNAME)
            NEW_DIRNAME = DIRNAME.replace(OLD_NAME, NEW_NAME)
            NEW_DIR_PATH = os.path.join(DIRPATH, NEW_DIRNAME)
            os.rename(OLD_DIR_PATH, NEW_DIR_PATH)
            print(f"Renamed directory: {OLD_DIR_PATH} -> {NEW_DIR_PATH}")

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="Recursively rename specific files, folders, and content containing the old_name.")
    PARSER.add_argument("DIRECTORY", help="The root directory to start renaming from.")
    PARSER.add_argument("OLD_NAME", help="The old string to replace.")
    PARSER.add_argument("NEW_NAME", help="The new string to replace with.")
    
    ARGS = PARSER.parse_args()

    RENAME_FILES_FOLDERS(ARGS.DIRECTORY, ARGS.OLD_NAME, ARGS.NEW_NAME)

