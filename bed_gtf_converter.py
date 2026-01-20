#!/usr/bin/env python3

import argparse
import os
import gzip
import zipfile

def read_file(filename):
    """Read a BED or GTF file, compressed (.gz/.zip) or uncompressed."""
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    elif filename.endswith('.zip'):
        with zipfile.ZipFile(filename, 'r') as z:
            first_file = z.namelist()[0]
            with z.open(first_file) as f:
                return [line.decode().strip() for line in f if line.decode().strip() and not line.decode().startswith('#')]
    else:
        with open(filename, 'r') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]

def write_file(lines, filename, compression=None):
    """Write lines to a file, optionally compressed as gz or zip."""
    if compression == 'gz':
        with gzip.open(filename, 'wt') as f:
            f.write('\n'.join(lines) + '\n')
    elif compression == 'zip':
        zipname = filename if filename.endswith('.zip') else filename + '.zip'
        inner_name = os.path.basename(filename).replace('.zip','')
        with zipfile.ZipFile(zipname, 'w', zipfile.ZIP_DEFLATED) as z:
            z.writestr(inner_name, '\n'.join(lines) + '\n')
    else:
        with open(filename, 'w') as f:
            f.write('\n'.join(lines) + '\n')

def bed_to_gtf(bed_lines):
    """Convert BED lines to GTF format, preserving multiple identifiers in the name field."""
    gtf_lines = []
    for line in bed_lines:
        fields = line.split('\t')
        if len(fields) < 3:
            continue
        chrom, start, end = fields[0], fields[1], fields[2]
        name = fields[3] if len(fields) > 3 else '.'
        score = fields[4] if len(fields) > 4 else '.'
        strand = fields[5] if len(fields) > 5 else '.'

        # Split the BED name field if it contains multiple identifiers
        gene_id = transcript_id = gene_name = name
        if '|' in name:
            parts = name.split('|')
            gene_id = parts[0]
            transcript_id = parts[1] if len(parts) > 1 else parts[0]
            gene_name = parts[2] if len(parts) > 2 else parts[0]

        gtf_fields = [
            chrom,
            "converted",
            "exon",
            str(int(start)+1),
            end,
            score,
            strand,
            ".",
            f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_name "{gene_name}";'
        ]
        gtf_lines.append('\t'.join(gtf_fields))
    return gtf_lines

def gtf_to_bed(gtf_lines):
    """Convert GTF lines to BED format, preserving gene_id, transcript_id, gene_name in the name field."""
    bed_lines = []
    for line in gtf_lines:
        fields = line.split('\t')
        if len(fields) < 9:
            continue
        chrom = fields[0]
        start = str(int(fields[3]) - 1)  # BED is 0-based
        end = fields[4]
        score = fields[5] if len(fields) > 5 else '.'
        strand = fields[6] if len(fields) > 6 else '.'
        attrs = fields[8]

        gene_id = transcript_id = gene_name = "."
        for attr in attrs.split(';'):
            attr = attr.strip()
            if attr.startswith('gene_id'):
                gene_id = attr.split(' ')[1].replace('"','')
            elif attr.startswith('transcript_id'):
                transcript_id = attr.split(' ')[1].replace('"','')
            elif attr.startswith('gene_name'):
                gene_name = attr.split(' ')[1].replace('"','')

        # Combine identifiers into BED name field
        name = '|'.join([gene_id, transcript_id, gene_name])
        bed_lines.append('\t'.join([chrom, start, end, name, score, strand]))
    return bed_lines

def main():
    parser = argparse.ArgumentParser(description="Convert between BED and GTF (supports gz/zip).")
    parser.add_argument('-i', '--input', required=True, help="Input BED or GTF file (can be .gz or .zip)")
    parser.add_argument('-o', '--output', required=True, help="Output file name (extension determines BED/GTF)")
    parser.add_argument('--compress', choices=['gz','zip','none'], default='none', help="Compression for output (gz, zip, none)")
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    compression = args.compress if args.compress != 'none' else None

    # Detect input format from filename
    if input_file.endswith('.bed') or '.bed.' in input_file:
        input_format = 'BED'
    elif input_file.endswith('.gtf') or '.gtf.' in input_file:
        input_format = 'GTF'
    else:
        raise ValueError("Cannot determine input format from filename. Use .bed or .gtf")

    # Detect output format from filename
    if output_file.endswith('.bed') or '.bed.' in output_file:
        output_format = 'BED'
    elif output_file.endswith('.gtf') or '.gtf.' in output_file:
        output_format = 'GTF'
    else:
        raise ValueError("Cannot determine output format from filename. Use .bed or .gtf")

    # Read input file
    lines = read_file(input_file)

    # Perform conversion
    if input_format == 'BED' and output_format == 'GTF':
        out_lines = bed_to_gtf(lines)
    elif input_format == 'GTF' and output_format == 'BED':
        out_lines = gtf_to_bed(lines)
    elif input_format == output_format:
        out_lines = lines  # same format, just handle compression
    else:
        raise ValueError("Unsupported conversion")

    # Write output file
    write_file(out_lines, output_file, compression)
    print(f"Converted {input_file} ({input_format}) -> {output_file} ({output_format}, compression={compression})")

if __name__ == "__main__":
    main()

# =============================================================================
# SCRIPT FUNCTIONALITY SUMMARY (ENHANCED)
# =============================================================================
# This Python utility script converts between BED and GTF formats with support for compressed files.
# Key features:
# 1. Automatically detects input format (BED or GTF) from the filename.
# 2. Supports input files compressed as .gz or .zip and reads them transparently.
# 3. Converts BED → GTF and GTF → BED:
#       - BED → GTF: Converts BED fields to GTF fields and splits the name field into gene_id, transcript_id, gene_name if '|' is used.
#       - GTF → BED: Combines gene_id, transcript_id, and gene_name into the BED name field using '|'.
# 4. Preserves key attributes to maintain identifiers across conversions.
# 5. Allows user to specify output compression: gz, zip, or none.
# 6. Writes output in requested format and compression, naming the inner file correctly if zipped.
# 7. Command-line interface provided with argparse.
# Example usage:
#   python bed_gtf_converter.py -i input.bed.gz -o output.gtf.gz --compress gz
#   python bed_gtf_converter.py -i input.gtf -o output.bed.zip --compress zip
# =============================================================================

