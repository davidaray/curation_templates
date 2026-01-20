import argparse
import csv

def parse_file1_line(line):
    parts = line.strip().split('\t')
    identifier, seq, coord, _, location = parts
    location_parts = location.split(':')
    chr_value = location_parts[1]
    coord_value = int(location_parts[2].split('-')[0]) - 15
    return chr_value, coord_value, seq

def load_file2(file2_path):
    file2_data = {}
    with open(file2_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            chr_value, coord = row[0], int(row[1])
            file2_data[(chr_value, coord)] = row
    return file2_data

def process_files(file1_path, file2_path, output_path):
    file2_data = load_file2(file2_path)
    output_rows = []
    
    with open(file1_path, 'r') as f1:
        for line in f1:
            chr_value, coord_value, seq = parse_file1_line(line)
            
            if (chr_value, coord_value) in file2_data:
                file2_row = file2_data[(chr_value, coord_value)]
                output_row = file2_row[:4] + [file2_row[5], "?", seq, "?", "?", "?", "?", "NONE", "INSIDER", "?", "?", "DELETION"]
                output_rows.append(output_row)
    
    with open(output_path, 'w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerows(output_rows)

def main():
    parser = argparse.ArgumentParser(description="Process TSD and deletion files.")
    parser.add_argument('-t', '--tsd', required=True, help="TSD input file")
    parser.add_argument('-d', '--deletion', required=True, help="Deletion input file")
    parser.add_argument('-o', '--output', required=True, help="Output TSV file")
    args = parser.parse_args()
    
    process_files(args.tsd, args.deletion, args.output)
    
if __name__ == "__main__":
    main()
