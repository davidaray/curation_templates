import argparse
import gzip
import subprocess

def filter_repeatmasker(input_file, min_div, max_div, output_file, threads):
    # Open input using gzip
    with gzip.open(input_file, 'rt') as infile:
        # Open pigz subprocess for parallel compression
        pigz_cmd = ['pigz', f'-p{threads}', '-c']
        with subprocess.Popen(pigz_cmd, stdin=subprocess.PIPE, stdout=open(output_file, 'wb'), universal_newlines=True) as pigz_proc:
            for line in infile:
                if line.startswith("#") or not line.strip():
                    continue  # Skip headers and empty lines

                fields = line.strip().split()
                if len(fields) < 2:
                    continue  # Skip malformed lines

                try:
                    divergence = float(fields[1])
                except ValueError:
                    continue  # Skip if divergence is not a float

                if min_div <= divergence <= max_div:
                    pigz_proc.stdin.write(line)
            pigz_proc.stdin.close()
            pigz_proc.wait()

def main():
    parser = argparse.ArgumentParser(description="Filter RepeatMasker .out.gz file by divergence.")
    parser.add_argument("-i", "--input", required=True, help="Input compressed RepeatMasker .out.gz file")
    parser.add_argument("-m", "--min_div", type=float, default=0.0, help="Minimum divergence threshold (default: 0.0)")
    parser.add_argument("-M", "--max_div", type=float, required=True, help="Maximum divergence threshold")
    parser.add_argument("-o", "--output", required=True, help="Output compressed .gz file for filtered results")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use for pigz compression")
    args = parser.parse_args()

    filter_repeatmasker(args.input, args.min_div, args.max_div, args.output, args.threads)

if __name__ == "__main__":
    main()


"""
Use python and argparse to:

1. take in a compressed RepeatMasker .out file (-i)
2. filter for maximum (-M) and minimum (-m) divergence thresholds
3. save the filtered file as a new file with a user-designated name (-o)

"""
