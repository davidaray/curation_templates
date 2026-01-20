import argparse
import gzip
import re
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey

def parse_repeatmasker_out(file):
    data = []
    with (gzip.open(file, 'rt') if file.endswith('.gz') else open(file, 'r')) as f:
        for line in f:
            if re.match(r'\s*\d+', line):  # Skip header and blank lines
                fields = line.split()
                repeat_class_family = fields[10]
                if '/' in repeat_class_family:
                    repeat_class, repeat_family = repeat_class_family.split('/')
                else:
                    repeat_class, repeat_family = repeat_class_family, 'Unknown'
                data.append({
                    'class': repeat_class,
                    'family': repeat_family
                })
    return data

def categorize_repeats(data):
    counts = defaultdict(float)
    for entry in data:
        repeat_class = entry['class']
        repeat_family = entry['family']

        # Classify the repeats
        if repeat_class == 'LINE':
            counts['Class I'] += 1
            counts['Non-LTR retrotransposon'] += 1
            counts[repeat_family] += 1
        elif repeat_class == 'LTR':
            counts['Class I'] += 1
            counts['LTR retrotransposon'] += 1
            counts[repeat_family] += 1
        elif repeat_class == 'SINE':
            counts['Class I'] += 1
            counts['Non-LTR retrotransposon'] += 1
        elif repeat_class == 'DNA':
            counts['Class II'] += 1
            counts['TIR TE'] += 1
        elif repeat_class == 'Helitron':
            counts['Class II'] += 1
            counts['Rolling circle'] += 1
        elif repeat_class == 'Simple_repeat':
            counts['Simple'] += 1
        elif repeat_class == 'Satellite':
            counts['Satellite'] += 1
        else:
            counts['Other'] += 1

    # Calculate Repetitive and Unmasked
    total_repeats = sum(counts.values())
    counts['Repetitive'] = total_repeats
    counts['Unmasked'] = 100 - total_repeats  # Assuming genome total is 100%
    
    return counts

def generate_sankey(counts):
    # Normalize counts to percentages relative to the total genome
    total_genome = counts['Repetitive'] + counts['Unmasked']
    
    # Scale all counts to sum up to 100
    for key in counts:
        counts[key] = (counts[key] / total_genome) * 100

    sankey = Sankey(unit=None)

    # Starting from Genome -> Repetitive and Unmasked
    sankey.add(flows=[counts['Repetitive'], counts['Unmasked']],
               labels=['Repetitive', 'Unmasked'], orientations=[1, -1])

    # Subdivide Repetitive into TEs and other categories
    sankey.add(flows=[counts['Transposable element'], counts['Satellite'], counts['Simple']],
               labels=['Transposable element', 'Satellite', 'Simple'],
               orientations=[1, -1, -1], prior=0, connect=(0, 0))

    # Subdivide TEs into Class I and Class II
    sankey.add(flows=[counts['Class I'], counts['Class II']],
               labels=['Class I', 'Class II'], orientations=[1, -1], prior=1, connect=(0, 0))

    # Subdivide Class I
    sankey.add(flows=[counts['Non-LTR retrotransposon'], counts['LTR retrotransposon'], counts['Penelope']],
               labels=['Non-LTR retrotransposon', 'LTR retrotransposon', 'Penelope'],
               orientations=[1, -1, -1], prior=2, connect=(0, 0))

    # Subdivide Class II
    sankey.add(flows=[counts['TIR TE'], counts['Rolling circle']],
               labels=['TIR TE', 'Rolling circle'], orientations=[1, -1], prior=3, connect=(0, 0))

    fig, ax = plt.subplots()
    sankey.finish()
    plt.title('Sankey Diagram of Repetitive Elements')
    plt.savefig(args.sankey_output)  # Save the output figure
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Generate a Sankey diagram from RepeatMasker output.')
    parser.add_argument('-r', '--repeatmasker', required=True, help='RepeatMasker .out file (can be gzipped)')
    parser.add_argument('-s', '--sankey', required=False, help='Output Sankey diagram')

    args = parser.parse_args()

    # Parse the RepeatMasker file
    repeat_data = parse_repeatmasker_out(args.repeatmasker)

    # Categorize the repeat elements
    counts = categorize_repeats(repeat_data)

    # Generate the Sankey diagram
    generate_sankey(counts)

if __name__ == "__main__":
    main()

