import sys
import random

def transform_clusters(input_file, output_file):
    sequences = []
    
    with open(input_file, 'r') as infile:
        current_cluster = None
        for line in infile:
            line = line.strip()
            if line.startswith("CLUSTER"):
                current_cluster = line.split()[1]
            elif line:
                sequences.append(f"{current_cluster} {line}")
    
    # Shuffle the sequences
    random.shuffle(sequences)
    
    with open(output_file, 'w') as outfile:
        for sequence in sequences:
            outfile.write(sequence + "\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python transform_clusters.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    transform_clusters(input_file, output_file)

if __name__ == '__main__':
    main()
