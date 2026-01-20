import argparse
import csv
from collections import defaultdict

def parse_tsv(input_file):
    te_class_counts = defaultdict(int)
    te_annotation_counts = defaultdict(lambda: [0, ""])
    
    all_te_classes = {"LINE", "SINE", "LTR", "DNA", "RC"}  # Ensure all are included
    total_polymorphisms = 0
    
    with open(input_file, "r") as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        
        for row in reader:
            total_polymorphisms += 1
            te_class = row["CLASS/FAMILY"].split("/")[0]  # Extract TE class
            te_class_counts[te_class] += 1
            
            te_annotation = row["TE_annotation"]
            te_annotation_counts[te_annotation][0] += 1
            te_annotation_counts[te_annotation][1] = row["CLASS/FAMILY"]
    
    return total_polymorphisms, te_class_counts, te_annotation_counts, all_te_classes

def write_output1(output_file, total_polymorphisms, te_class_counts, all_te_classes):
    with open(output_file, "w") as out:
        out.write("TE_Class\tCount\n")
        out.write(f"Total\t{total_polymorphisms}\n")
        
        for te_class in sorted(all_te_classes):
            out.write(f"{te_class}\t{te_class_counts.get(te_class, 0)}\n")

def write_output2(output_file, te_annotation_counts):
    sorted_annotations = sorted(te_annotation_counts.items(), key=lambda x: x[1][0], reverse=True)
    
    with open(output_file, "w") as out:
        out.write("TE_annotation\tCLASS/FAMILY\tCount\n")
        for annotation, (count, class_family) in sorted_annotations:
            out.write(f"{annotation}\t{class_family}\t{count}\n")

def main():
    parser = argparse.ArgumentParser(description="Process polymorphism data")
    parser.add_argument("-i", required=True, help="Input .tsv file")
    parser.add_argument("-o1", required=True, help="Output file for TE class counts")
    parser.add_argument("-o2", required=True, help="Output file for TE annotation counts")
    
    args = parser.parse_args()
    
    total_polymorphisms, te_class_counts, te_annotation_counts, all_te_classes = parse_tsv(args.i)
    write_output1(args.o1, total_polymorphisms, te_class_counts, all_te_classes)
    write_output2(args.o2, te_annotation_counts)

if __name__ == "__main__":
    main()
