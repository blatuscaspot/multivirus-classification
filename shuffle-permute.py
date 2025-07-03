import random

def read_fastq(file_path):
    """Yield tuples of (header, seq, plus, qual)."""
    with open(file_path, 'r') as f:
        while True:
            lines = [f.readline().rstrip() for _ in range(4)]
            if not lines[0]:
                break
            yield tuple(lines)

def write_fastq(records, output_path):
    with open(output_path, 'w') as out_f:
        for record in records:
            out_f.write('\n'.join(record) + '\n')

def shuffle_reads(records):
    """Return list of records with shuffled sequence and quality lines (within each read)."""
    shuffled = []
    for h, seq, p, qual in records:
        seq_list = list(seq)
        qual_list = list(qual)
        zipped = list(zip(seq_list, qual_list))
        random.shuffle(zipped)
        shuffled_seq, shuffled_qual = zip(*zipped)
        shuffled.append((h, ''.join(shuffled_seq), p, ''.join(shuffled_qual)))
    return shuffled

def permute_headers(records):
    """Return list of records with headers randomly reassigned to other reads."""
    headers = [r[0] for r in records]
    random.shuffle(headers)
    permuted = [(new_h, r[1], r[2], r[3]) for new_h, r in zip(headers, records)]
    return permuted

# List of files
input_files = ["Original.fastq", "Enriched.fastq", "Depleted.fastq"]

for input_file in input_files:
    print(f"Processing {input_file}...")
    records = list(read_fastq(input_file))

    # Shuffled version
    shuffled_records = shuffle_reads(records)
    shuffled_output = input_file.replace('.fastq', '_shuffled.fastq')
    write_fastq(shuffled_records, shuffled_output)
    print(f"→ Wrote shuffled file: {shuffled_output}")

    # Permuted version
    permuted_records = permute_headers(records)
    permuted_output = input_file.replace('.fastq', '_permuted.fastq')
    write_fastq(permuted_records, permuted_output)
    print(f"→ Wrote permuted file: {permuted_output}")

print("\nDone!")
