#!/usr/bin/env python3
"""
Resolve overlapping regions in RepeatMasker output files.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict


class RepeatElement:
    """Represents a single RepeatMasker annotation."""
    
    def __init__(self, line, line_num):
        self.original_line = line
        self.line_num = line_num
        parts = line.rstrip('\n').split()
        
        # Check for trailing *
        self.has_star = parts[-1] == '*'
        if self.has_star:
            parts = parts[:-1]
        
        # Parse columns
        self.score = int(parts[0])
        self.divergence = float(parts[1])
        self.deletion = float(parts[2])
        self.insertion = float(parts[3])
        self.sequence = parts[4]
        self.begin = int(parts[5])
        self.end = int(parts[6])
        self.left = parts[7]
        self.strand = parts[8]
        self.repeat_name = parts[9]
        self.repeat_class = parts[10]
        self.repeat_begin = parts[11]
        self.repeat_end = parts[12]
        self.repeat_left = parts[13]
        self.id = int(parts[14])
        
        # Check if it's a Simple Repeat
        self.is_simple_repeat = 'Simple_repeat' in self.repeat_class or 'Low_complexity' in self.repeat_class
        
    def length(self):
        return self.end - self.begin + 1
    
    def overlaps_with(self, other):
        """Check if this element overlaps with another."""
        if self.sequence != other.sequence:
            return False
        return not (self.end < other.begin or self.begin > other.end)
    
    def contains(self, other):
        """Check if this element completely contains another."""
        if self.sequence != other.sequence:
            return False
        return self.begin <= other.begin and self.end >= other.end
    
    def get_effective_divergence(self):
        """Get divergence, treating simple repeats as highest."""
        if self.is_simple_repeat:
            return 999.9  # Very high divergence
        return self.divergence


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Resolve overlapping regions in RepeatMasker output'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input RepeatMasker .out file'
    )
    parser.add_argument(
        '-s', '--split',
        choices=['class', 'family', 'none'],
        default='none',
        help='Split output by TE Class or Family (default: none)'
    )
    parser.add_argument(
        '-ot', '--output-type',
        required=True,
        choices=['bed', 'out', 'both'],
        help='Output file type'
    )
    parser.add_argument(
        '-p', '--prefix',
        required=True,
        help='Output file prefix'
    )
    parser.add_argument(
        '-ov', '--overlap-resolution',
        choices=['higher_score', 'longer_element', 'lower_divergence'],
        default='lower_divergence',
        help='Method for resolving overlaps (default: lower_divergence)'
    )
    parser.add_argument(
        '-m', '--min-hits',
        type=int,
        default=0,
        help='Minimum number of hits per TE Class/Family to include (default: 0)'
    )
    parser.add_argument(
        '-dmax', '--max-divergence',
        type=float,
        help='Maximum divergence allowed'
    )
    parser.add_argument(
        '-dmin', '--min-divergence',
        type=float,
        help='Minimum divergence allowed'
    )
    return parser.parse_args()


def read_repeatmasker_file(filepath):
    """Read RepeatMasker file and return header and elements."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Store header (first 3 lines)
    header = ''.join(lines[:3])
    
    # Parse elements
    elements = []
    for i, line in enumerate(lines[3:], start=4):
        if line.strip():
            try:
                elements.append(RepeatElement(line, i))
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line {i}: {line.strip()}", file=sys.stderr)
    
    return header, elements


def filter_by_divergence(elements, min_div, max_div):
    """Filter elements by divergence range."""
    filtered = []
    for elem in elements:
        # Don't filter simple repeats by divergence
        if elem.is_simple_repeat:
            filtered.append(elem)
            continue
        
        if min_div is not None and elem.divergence < min_div:
            continue
        if max_div is not None and elem.divergence > max_div:
            continue
        filtered.append(elem)
    
    return filtered


def resolve_overlap(elem1, elem2, method):
    """
    Determine which element to keep in case of overlap.
    Returns the element to keep, or None if they should be split.
    """
    # Check for containment
    elem1_contains_elem2 = elem1.contains(elem2)
    elem2_contains_elem1 = elem2.contains(elem1)
    
    if elem1_contains_elem2 or elem2_contains_elem1:
        # Containment case
        if method == 'higher_score':
            if elem1.score > elem2.score:
                return elem1
            elif elem2.score > elem1.score:
                return elem2
            else:
                return elem1  # First if equal
        
        elif method == 'longer_element':
            if elem1.length() > elem2.length():
                return elem1
            elif elem2.length() > elem1.length():
                return elem2
            else:
                return elem1  # First if equal
        
        elif method == 'lower_divergence':
            div1 = elem1.get_effective_divergence()
            div2 = elem2.get_effective_divergence()
            if div1 < div2:
                return elem1
            elif div2 < div1:
                return elem2
            else:
                return elem1  # First if equal
    
    else:
        # Basic overlap case - return element to prioritize
        if method == 'higher_score':
            return elem1 if elem1.score >= elem2.score else elem2
        elif method == 'longer_element':
            return elem1 if elem1.length() >= elem2.length() else elem2
        elif method == 'lower_divergence':
            div1 = elem1.get_effective_divergence()
            div2 = elem2.get_effective_divergence()
            return elem1 if div1 <= div2 else elem2
    
    return elem1


def resolve_overlaps(elements, method):
    """Resolve all overlaps in the element list."""
    # Group by sequence
    by_sequence = defaultdict(list)
    for elem in elements:
        by_sequence[elem.sequence].append(elem)
    
    resolved = []
    
    for sequence, seq_elements in by_sequence.items():
        # Sort by start position
        seq_elements.sort(key=lambda e: (e.begin, e.end))
        
        # Process overlaps
        kept = []
        for elem in seq_elements:
            # Check if this element overlaps with any kept element
            overlaps = []
            for i, kept_elem in enumerate(kept):
                if elem.overlaps_with(kept_elem):
                    overlaps.append(i)
            
            if not overlaps:
                # No overlap, keep it
                kept.append(elem)
            else:
                # Has overlaps - resolve them
                should_keep = True
                indices_to_remove = []
                
                for idx in overlaps:
                    kept_elem = kept[idx]
                    winner = resolve_overlap(kept_elem, elem, method)
                    
                    # Check for containment
                    if kept_elem.contains(elem) or elem.contains(kept_elem):
                        # Containment - keep only the winner
                        if winner == elem:
                            indices_to_remove.append(idx)
                        else:
                            should_keep = False
                            break
                    else:
                        # Basic overlap - trim the loser
                        if winner == elem:
                            # Trim the kept element to not overlap
                            indices_to_remove.append(idx)
                        else:
                            # Don't add current element
                            should_keep = False
                            break
                
                # Remove elements that lost
                for idx in sorted(indices_to_remove, reverse=True):
                    kept.pop(idx)
                
                if should_keep:
                    kept.append(elem)
        
        resolved.extend(kept)
    
    return resolved


def write_out_format(elements, filepath, header):
    """Write elements in RepeatMasker .out format."""
    with open(filepath, 'w') as f:
        f.write(header)
        for elem in elements:
            f.write(elem.original_line)


def write_bed_format(elements, filepath):
    """Write elements in BED format."""
    with open(filepath, 'w') as f:
        for elem in elements:
            # BED format: chrom, start, end, name, score, strand
            strand = '+' if elem.strand == '+' else '-'
            name = f"{elem.repeat_name}#{elem.repeat_class}"
            f.write(f"{elem.sequence}\t{elem.begin-1}\t{elem.end}\t{name}\t{elem.score}\t{strand}\n")


def main():
    """Main processing pipeline."""
    args = parse_arguments()
    
    print(f"Reading {args.input}...")
    header, elements = read_repeatmasker_file(args.input)
    print(f"Read {len(elements)} elements")
    
    # Filter by divergence if specified
    if args.min_divergence is not None or args.max_divergence is not None:
        elements = filter_by_divergence(elements, args.min_divergence, args.max_divergence)
        print(f"After divergence filtering: {len(elements)} elements")
    
    # Resolve overlaps
    print(f"Resolving overlaps using '{args.overlap_resolution}' method...")
    resolved = resolve_overlaps(elements, args.overlap_resolution)
    print(f"After overlap resolution: {len(resolved)} elements")
    
    # Apply minimum hits filter if splitting
    if args.split != 'none':
        # Count hits per class/family
        counts = defaultdict(int)
        for elem in resolved:
            key = elem.repeat_class if args.split == 'class' else elem.repeat_class.split('/')[0]
            counts[key] += 1
        
        # Filter
        filtered = []
        for elem in resolved:
            key = elem.repeat_class if args.split == 'class' else elem.repeat_class.split('/')[0]
            if counts[key] >= args.min_hits:
                filtered.append(elem)
        resolved = filtered
        print(f"After minimum hits filter: {len(resolved)} elements")
    
    # Output files
    if args.split == 'none':
        # Single output file
        if args.output_type in ['out', 'both']:
            out_file = f"{args.prefix}.out"
            write_out_format(resolved, out_file, header)
            print(f"Wrote {out_file}")
        
        if args.output_type in ['bed', 'both']:
            bed_file = f"{args.prefix}.bed"
            write_bed_format(resolved, bed_file)
            print(f"Wrote {bed_file}")
    
    else:
        # Split by class or family
        by_category = defaultdict(list)
        for elem in resolved:
            if args.split == 'class':
                category = elem.repeat_class
            else:  # family
                category = elem.repeat_class.split('/')[0] if '/' in elem.repeat_class else elem.repeat_class
            by_category[category].append(elem)
        
        print(f"\nWriting {len(by_category)} split files...")
        for category, cat_elements in by_category.items():
            safe_cat = category.replace('/', '_')
            
            if args.output_type in ['out', 'both']:
                out_file = f"{args.prefix}.{safe_cat}.out"
                write_out_format(cat_elements, out_file, header)
                print(f"  {safe_cat}: {len(cat_elements)} elements -> {out_file}")
            
            if args.output_type in ['bed', 'both']:
                bed_file = f"{args.prefix}.{safe_cat}.bed"
                write_bed_format(cat_elements, bed_file)
                print(f"  {safe_cat}: {len(cat_elements)} elements -> {bed_file}")
    
    print("\nDone!")


if __name__ == '__main__':
    main()
