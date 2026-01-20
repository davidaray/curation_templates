import os
import argparse
import shutil
import fnmatch

def delete_files_and_folders(root_dir):
    # Define the files and folders to be deleted
    folder_patterns = ['assemblies_dir', 'rmasker_dir', 'RMPart', 'RM_*']
    file_patterns = ['*.fa.gz', '*.log', '*.nhr', '*.nin', '*.njs', '*.nnd', '*.nni', '*.nog', '*.nsq', '*.RMrun.out', '*.translation']
    exclude_files = '*.softmasked.fa.gz'

    # Walk through the directory
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # Delete matching folders
        for folder_pattern in folder_patterns:
            for dirname in fnmatch.filter(dirnames, folder_pattern):
                folder_path = os.path.join(dirpath, dirname)
                print(f"Deleting folder: {folder_path}")
                shutil.rmtree(folder_path)

        # Delete matching files
        for file_pattern in file_patterns:
            for filename in fnmatch.filter(filenames, file_pattern):
                if not fnmatch.fnmatch(filename, exclude_files):
                    file_path = os.path.join(dirpath, filename)
                    print(f"Deleting file: {file_path}")
                    os.remove(file_path)

if __name__ == '__main__':
    # Set up argparse to get the directory from the command line
    parser = argparse.ArgumentParser(description="Delete specific files and folders.")
    parser.add_argument('directory', type=str, help='Directory to search and delete files/folders')
    
    args = parser.parse_args()

    # Call the function with the provided directory
    delete_files_and_folders(args.directory)

