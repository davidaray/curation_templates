#!/usr/bin/env python3

import argparse
import gzip
import re

def open_file(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def parse_repeatmasker_line(line):
    # RepeatMasker lines have a fixed-width format. Use regex to parse.
    fields = re.split(r'\s+', line.strip())

    # Example parsed fields:
    # 0  - SW score
    # 1  - % div.
    # 2  - % del.
    # 3  - % ins.
    # 4  - query sequence
    # 5  - query start
    # 6  - query end
    # 7  - query remaining
    # 8  - strand
    # 9  - matching repeat
    # 10 - repeat class/family
    # 11 - repeat start
    # 12 - repeat end
    # 13 - repeat left
    # 14 - ID

    if len(fields) < 14:
        return None

    sw_score = fields[0]
    chrom = fields[4]
    start = int(fields[5])
    end = int(fields[6])
    strand = '-' if fields[8] == 'C' else '+'
    repeat_name = fields[9]
    repeat_class = fields[10]
    repeat_start = fields[11]
    repeat_end = fields[12]
    repeat_left = fields[13]

    attributes = f"ID={repeat_name}_{start}_{end};Name={repeat_name};Class={repeat_class};SW_score={sw_score}"

    return {
        'seqid': chrom,
        'source': 'RepeatMasker',
        'type': 'repeat_region',
        'start': start,
        'end': end,
        'score': sw_score,
        'strand': strand,
        'phase': '.',
        'attributes': attributes
    }

def convert_repeatmasker_to_gff(input_file, output_file):
    with open_file(input_file) as infile, open(output_file, 'w') as out:
        out.write("##gff-version 3\n")
        started = False
        for line in infile:
            if not started:
                if line.startswith("   SW"):
                    started = True
                continue
            if line.strip() == "" or line.startswith("score"):
                continue

            parsed = parse_repeatmasker_line(line)
            if parsed:
                out.write("\t".join([
                    parsed['seqid'],
                    parsed['source'],
                    parsed['type'],
                    str(parsed['start']),
                    str(parsed['end']),
                    parsed['score'],
                    parsed['strand'],
                    parsed['phase'],
                    parsed['attributes']
                ]) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Convert RepeatMasker .out or .out.gz file to GFF3 format.")
    parser.add_argument("-i", "--input", required=True, help="RepeatMasker .out or .out.gz file")
    parser.add_argument("-o", "--output", required=True, help="Output GFF3 file")
    args = parser.parse_args()

    convert_repeatmasker_to_gff(args.input, args.output)

if __name__ == "__main__":
    main()

