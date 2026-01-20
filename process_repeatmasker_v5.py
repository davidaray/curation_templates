#!/usr/bin/env python3
"""
resolve_rm_overlaps.py

Resolve overlapping regions in RepeatMasker .out files using user-specified rules.

Behavior summary (Option A - trimming):
 - Basic overlap (partial): overlapping bases are assigned to the winner; the losing
   element is trimmed to remove the overlapping bases (may be split into two fragments).
 - Containment: the losing element (by the chosen rule) is discarded entirely.
 - Simple repeats: treated as having a very large divergence so they lose in
   'lower_divergence' mode.

Author: ChatGPT
"""

from __future__ import annotations
import argparse
import re
import sys
from collections import defaultdict, namedtuple
from dataclasses import dataclass, field
from typing import List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor

print("DEBUG: script started")

# Data structure for an interval / element
@dataclass
class Element:
    idx: int                      # original line index (for stable tie-break)
    raw_line: str                 # original full text line
    score: float
    div: float                    # divergence (float); Simple repeats get large value
    seqname: str                  # query sequence (chrom/scaffold)
    start: int                    # begin
    end: int                      # end (inclusive as in RepeatMasker)
    strand: str
    repeat_name: str
    class_family: str             # e.g., "SINE/Ves" or "LTR/ERV"
    extra_tokens: List[str] = field(default_factory=list)  # remaining tokens after class/fam
    source_order: int = 0         # original order for tie-breaking

    def length(self) -> int:
        return self.end - self.start + 1

    def copy_with_coords(self, new_start: int, new_end: int, new_idx: Optional[int] = None):
        return Element(
            idx=(new_idx if new_idx is not None else self.idx),
            raw_line=self.raw_line,
            score=self.score,
            div=self.div,
            seqname=self.seqname,
            start=new_start,
            end=new_end,
            strand=self.strand,
            repeat_name=self.repeat_name,
            class_family=self.class_family,
            extra_tokens=list(self.extra_tokens),
            source_order=self.source_order
        )

def parse_rm_line(line: str, idx: int, simple_divergence_value: float = 1e6) -> Optional[Element]:
    """
    Parse a single RepeatMasker .out line according to the header layout you provided:
       SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
    score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

    We parse tokens, locate the class/family token (it contains '/'), and reconstruct fields.
    Returns Element or None (for header/blank lines).
    """
    raw = line.rstrip("\n")
    if not raw.strip():
        return None
    # skip lines that look like header lines
    if raw.lstrip().startswith("SW") or raw.lstrip().startswith("score"):
        return None

    tokens = re.split(r'\s+', raw.strip())
    if len(tokens) < 10:
        # can't parse
        raise ValueError(f"Unable to parse line {idx+1}: '{raw}'")

    # standard mapping by token indices (empirical from your sample)
    # 0: score, 1: div, 2: del, 3: ins, 4: sequence, 5: begin, 6: end, 7: (left), 8: strand, 9..: repeat name & class/family ...
    try:
        score = float(tokens[0])
        div = tokens[1]
        # some divergence entries like "."? attempt float
        try:
            divf = float(div)
        except:
            divf = float('inf')
        seqname = tokens[4]
        start = int(tokens[5])
        end = int(tokens[6])
        # token 7 often like '(213170481)' -> ignore
        # strand:
        strand = tokens[8]
    except Exception as e:
        raise ValueError(f"Error parsing core columns on line {idx+1}: {e}\nLine: {raw}")

    # locate first token that includes a slash '/' which indicates class/family
    slash_idx = None
    for i in range(9, len(tokens)):
        if '/' in tokens[i]:
            slash_idx = i
            break
    if slash_idx is None:
        # fall back: last token that looks like a class (no spaces is ideal)
        # try last-4 token heuristics
        slash_idx = len(tokens) - 4 if len(tokens) >= 12 else len(tokens) - 2

    repeat_name = ' '.join(tokens[9:slash_idx])
    class_family = tokens[slash_idx]
    extra = tokens[slash_idx+1:]

    # Treat Simple repeats as having huge divergence so they lose in lower_divergence mode
    if 'simple' in class_family.lower() or 'satellite' in class_family.lower():
        divf = simple_divergence_value

    e = Element(
        idx=idx,
        raw_line=raw,
        score=score,
        div=divf,
        seqname=seqname,
        start=start,
        end=end,
        strand=strand,
        repeat_name=repeat_name,
        class_family=class_family,
        extra_tokens=extra,
        source_order=idx
    )
    return e

def overlaps(a: Element, b: Element) -> bool:
    return not (a.end < b.start or b.end < a.start)

def is_containment(a: Element, b: Element) -> bool:
    """Return True if a contains b or b contains a."""
    return (a.start <= b.start and a.end >= b.end) or (b.start <= a.start and b.end >= a.end)

def containment_holder(a: Element, b: Element) -> Optional[Tuple[Element, Element]]:
    """
    Return (container, containee) if one contains the other; else None.
    """
    if a.start <= b.start and a.end >= b.end:
        return (a, b)
    if b.start <= a.start and b.end >= a.end:
        return (b, a)
    return None

def choose_winner(a: Element, b: Element, mode: str) -> Element:
    """
    Choose winner element for overlapping region according to mode.
    Tie-breaking: smaller source_order (earlier in file) wins.
    """
    if mode == "higher_score":
        if a.score > b.score:
            return a
        if b.score > a.score:
            return b
    elif mode == "longer_element":
        if a.length() > b.length():
            return a
        if b.length() > a.length():
            return b
    elif mode == "lower_divergence":
        # lower numeric divergence wins
        if a.div < b.div:
            return a
        if b.div < a.div:
            return b
    # tie -> earlier in file (lower source_order)
    return a if a.source_order <= b.source_order else b

def trim_element_remove_overlap(loser: Element, overlap_start: int, overlap_end: int) -> List[Element]:
    """
    Trim the loser by removing the [overlap_start, overlap_end] segment.
    Returns list of 0, 1, or 2 new Elements corresponding to remaining fragments.
    """
    new_elements = []
    # left fragment
    if loser.start < overlap_start:
        left = loser.copy_with_coords(loser.start, overlap_start - 1)
        new_elements.append(left)
    # right fragment
    if loser.end > overlap_end:
        right = loser.copy_with_coords(overlap_end + 1, loser.end)
        new_elements.append(right)
    return new_elements

def process_intervals(elements: List[Element], mode: str = "lower_divergence") -> List[Element]:
    """
    Main function to resolve overlaps. Processes elements in original order and builds
    a list of resolved elements (may be trimmed or split).
    """
    resolved: List[Element] = []

    for e in elements:
        # We'll maintain a working queue of fragments of e to insert
        fragments: List[Element] = [e]
        i = 0
        while i < len(fragments):
            cur = fragments[i]
            j = 0
            while j < len(resolved):
                existing = resolved[j]
                if cur.seqname != existing.seqname:
                    j += 1
                    continue
                if not overlaps(cur, existing):
                    j += 1
                    continue

                # If containment, special rules: the lower-scoring (per mode) element is discarded entirely
                cont = containment_holder(cur, existing)
                if cont:
                    container, containee = cont
                    # pick winner between container and containee using mode (BUT tie-break rules as specified)
                    winner = choose_winner(container, containee, mode)
                    loser = containee if winner is container else container

                    # If loser is container (i.e., winner was containee), special handling:
                    # Our rule: discard the lower scoring element entirely (even if it's the containing alignment).
                    # So we remove loser entirely.
                    if loser is existing:
                        # remove existing from resolved
                        resolved.pop(j)
                        # do not increment j (we removed it)
                        continue
                    else:
                        # loser is current fragment -> discard and break out of checking existing overlaps
                        fragments.pop(i)
                        i -= 1
                        break
                else:
                    # basic partial overlap -> determine winner for overlapping bases
                    winner = choose_winner(cur, existing, mode)
                    loser = existing if winner is cur else cur
                    # compute overlap coordinates
                    ov_s = max(cur.start, existing.start)
                    ov_e = min(cur.end, existing.end)

                    if loser is existing:
                        # trim existing in resolved: replace existing by fragments (0/1/2)
                        new_frags = trim_element_remove_overlap(existing, ov_s, ov_e)
                        # replace resolved[j] by the new_frags (maintain stable order)
                        resolved.pop(j)
                        # insert new_frags at position j (in order)
                        for f in reversed(new_frags):
                            resolved.insert(j, f)
                        # since we modified resolved, continue checking same j with same cur
                        continue
                    else:
                        # loser is current fragment: trim cur into up to two fragments
                        new_fragments = trim_element_remove_overlap(cur, ov_s, ov_e)
                        # replace cur in fragments list by new_fragments
                        fragments.pop(i)
                        for f in reversed(new_fragments):
                            fragments.insert(i, f)
                        i -= 1  # because we will increment at loop bottom
                        break  # break to re-evaluate newly-inserted fragment at same index
                j += 1
            i += 1
        # after resolving against all existing resolved intervals, append remaining fragments to resolved
        # preserve original relative order by simply extending
        resolved.extend(fragments)

    # final pass: sort resolved by seqname then start then original source order
    resolved.sort(key=lambda x: (x.seqname, x.start, x.source_order))
    return resolved

# Output helpers
def element_to_bed(e: Element) -> str:
    # BED is 0-based, half-open: start0 = start-1, end = end
    start0 = e.start - 1
    name = f"{e.repeat_name}|{e.class_family}|idx{e.idx}"
    score = int(e.score) if e.score is not None else 0
    return f"{e.seqname}\t{start0}\t{e.end}\t{name}\t{score}\t{e.strand}"

def element_to_out_line(e: Element) -> str:
    """
    Reconstruct a RepeatMasker-like .out line. This is approximate: many columns such as
    (left) values aren't recalculated exactly. We preserve original tokens where possible.
    """
    # We'll attempt to use the original raw_line but replace the query begin/end and the repeat-position begin/end if present.
    # Simpler: fabricate a consistent line:
    # score div del ins sequence begin end (left) strand repeat class/family <extra...>
    left_placeholder = "(0)"
    tokens = [
        f"{int(e.score):6d}",
        f"{e.div:5.1f}" if e.div < 1e5 else "  .",  # if simple repeat had huge div, leave as '.'
        "  0.0", "  0.0",
        f"{e.seqname}",
        f"{e.start:6d}",
        f"{e.end:6d}",
        f"{left_placeholder}",
        f"{e.strand}",
        f"{e.repeat_name}",
        f"{e.class_family}"
    ]
    # append extras (like repeat begin/end etc) if known
    if e.extra_tokens:
        tokens.extend(e.extra_tokens)
    return " ".join(tokens)

def write_outputs(resolved: List[Element], prefix: str, out_types: List[str], split_by: Optional[str],
                  progress_dir: Optional[str], min_hits: int):
    """
    Write outputs. If split_by is 'class' or 'family', create per-group outputs as well
    and apply min_hits filtering per group.
    """
    # group elements
    groups = defaultdict(list)
    if split_by == "class":
        for e in resolved:
            group_key = e.class_family  # class/family token like "SINE/Ves" - use whole string
            groups[group_key].append(e)
    elif split_by == "family":
        # family is the left side of class/family token if it exists like "SINE/Ves"
        for e in resolved:
            fam = e.class_family.split('/')[0] if '/' in e.class_family else e.class_family
            groups[fam].append(e)
    else:
        groups["all"] = resolved

    # optionally create progress dir
    if progress_dir:
        import os
        os.makedirs(progress_dir, exist_ok=True)

    def write_group(key_and_list):
        key, elems = key_and_list
        # apply min_hits if requested
        if min_hits is not None and len(elems) < min_hits:
            return f"skipped:{key}"

        # filename prefix
        safe_key = key.replace('/', '_').replace(' ', '_')
        base = f"{prefix}.{safe_key}" if key != "all" else prefix

        outputs = []
        if "bed" in out_types:
            bed_path = f"{base}.bed"
            with open(bed_path, "w") as h:
                for e in elems:
                    h.write(element_to_bed(e) + "\n")
            outputs.append(bed_path)
        if "out" in out_types:
            out_path = f"{base}.out"
            with open(out_path, "w") as h:
                # write a header line similar to RepeatMasker .out
                h.write(" SW  perc perc perc  query  position in query   matching   repeat  position in repeat\n")
                for e in elems:
                    h.write(element_to_out_line(e) + "\n")
            outputs.append(out_path)

        if progress_dir:
            # copy to progress dir as a status indicator
            import shutil, os
            for p in outputs:
                shutil.copy2(p, os.path.join(progress_dir, os.path.basename(p)))
        return f"written:{key}:{len(elems)}"

    # parallelize writes if many groups
    items = list(groups.items())
    results = []
    if len(items) > 1:
        with ThreadPoolExecutor() as ex:
            for res in ex.map(write_group, items):
                results.append(res)
    else:
        for res in map(write_group, items):
            results.append(res)

    return results

def compute_group_counts(elements: List[Element]) -> Tuple[dict, dict]:
    class_counts = defaultdict(int)
    family_counts = defaultdict(int)
    for e in elements:
        class_counts[e.class_family] += 1
        fam = e.class_family.split('/')[0] if '/' in e.class_family else e.class_family
        family_counts[fam] += 1
    return class_counts, family_counts

def main():
    parser = argparse.ArgumentParser(description="Resolve overlaps in RepeatMasker .out files (trim losing bases).")
    parser.add_argument('-i', '--input', required=True, help='Input RepeatMasker .out file')
    parser.add_argument('-s', '--split', choices=['class', 'family'], default=None,
                        help='Split outputs by TE class or family (optional)')
    parser.add_argument('-ot', '--output-type', required=True, choices=['bed', 'out', 'both'],
                        help='Output type: bed, out, or both')
    parser.add_argument('-p', '--prefix', required=True, help='Output filename prefix')
    parser.add_argument('-ov', '--ovlp_resolution', choices=['higher_score', 'longer_element', 'lower_divergence'],
                        default='lower_divergence', help='Overlap resolution mode')
    parser.add_argument('-m', '--min_hits', type=int, default=None,
                        help='Minimum number of hits from a given TE Class/Family to include that group in output (applies when splitting)')
    parser.add_argument('-dmax', '--dmax', type=float, default=None, help='Maximum divergence allowed (filter)')
    parser.add_argument('-dmin', '--dmin', type=float, default=None, help='Minimum divergence allowed (filter)')
    parser.add_argument('--progress-dir', default=None, help='Optional progress directory to copy outputs into')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads (used for parallel writing)')
    args = parser.parse_args()

    print("DEBUG: starting multiprocessing with", args.threads, "threads")

    # read file and parse lines
    elements: List[Element] = []
    with open(args.input, 'r') as fh:
        for idx, line in enumerate(fh):
            try:
                el = parse_rm_line(line, idx)
            except ValueError:
                el = None
            if el:
                elements.append(el)

    if not elements:
        print("No valid RepeatMasker entries parsed. Exiting.", file=sys.stderr)
        sys.exit(1)

    # apply divergence filters if provided (note: Simple repeats have huge divergence)
    filtered_elements = []
    for e in elements:
        keep = True
        if args.dmin is not None and (e.div < args.dmin):
            keep = False
        if args.dmax is not None and (e.div > args.dmax):
            keep = False
        if keep:
            filtered_elements.append(e)
    elements = filtered_elements

    if not elements:
        print("No entries left after divergence filtering. Exiting.", file=sys.stderr)
        sys.exit(1)

    # if splitting and min_hits provided, we'll apply min_hits only after group counts computed on original dataset
    class_counts, family_counts = compute_group_counts(elements)
    if args.split == 'class' and args.min_hits is not None:
        # drop any elements from classes with too few hits
        keep_classes = {k for k,v in class_counts.items() if v >= args.min_hits}
        elements = [e for e in elements if e.class_family in keep_classes]
    if args.split == 'family' and args.min_hits is not None:
        keep_fams = {k for k,v in family_counts.items() if v >= args.min_hits}
        elements = [e for e in elements if (e.class_family.split('/')[0] if '/' in e.class_family else e.class_family) in keep_fams]

    # Process overlaps
    resolved = process_intervals(elements, mode=args.ovlp_resolution)

    # Output
    if args.output_type == 'both':
        out_types = ['bed','out']
    else:
        out_types = [args.output_type]

    results = write_outputs(resolved, args.prefix, out_types, args.split, args.progress_dir, args.min_hits)
    print("Write summary:")
    for r in results:
        print("  ", r)

if __name__ == "__main__":
    main()
