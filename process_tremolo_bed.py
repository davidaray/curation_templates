import argparse

def process_tremolo_bed(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        line_count = 0
        for line in infile:
            line = line.strip()
            line_count += 1

            # Skip the second line
            if line_count == 2:
                continue

            # Skip empty lines
            if not line:
                continue

            # Identify and write the header line
            if line.startswith("#chrom"):
                outfile.write("#chrom\tstart\tend\tTE_annotation\tAssemblytics_value\tstrand\tTSD\tpident\tpsize_TE\tSIZE_TE\tNEW_POS\tFREQ\tFREQ_WITH_CLIPPED\tSV_SIZE\tID_TrEMOLO\tTYPE\n")
                continue

            columns = line.split("\t")

            # Ensure minimum column count
            if len(columns) < 15:
                continue

            # Extract annotation and assemblytics value
            annotation_parts = columns[3].split('|')
            te_annotation = annotation_parts[0]
            assemblytics_value = annotation_parts[1] if len(annotation_parts) > 1 else ""

            # Construct the new line with correct column order
            new_line = columns[:3] + [te_annotation, assemblytics_value] + columns[4:]

            # Apply filters after columns are finalized
            tsd = new_line[6]  # TSD column
            id_tremolo = new_line[14]  # ID_TrEMOLO column

            # Skip lines based on TSD or ID_TrEMOLO
            if tsd == "NONE" or ".Repeat_expansion" in id_tremolo:
                continue

            # Write the modified line
            outfile.write("\t".join(new_line) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Process Tremolo BED file")
    parser.add_argument("-i", "--input", required=True, help="Input BED file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    process_tremolo_bed(args.input, args.output)


if __name__ == "__main__":
    main()
