import sys
import re
import math
from collections import defaultdict
from Bio import SeqIO

def load_sequences(fasta_file):
	db = {}
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            db[record.id] = record.seq
    return db

def process_file(f1, readlength, offset):
    final_data = []
    for line in f1:
        line = line.strip()
        if line.startswith("@"):
            continue
        array = line.split()
        read_len = [int(x) for x in readlength.split(',')]
        off_set = [int(x) for x in offset.split(',')]
        final_data.append((array, read_len, off_set))
    return final_data

def process_duplicates_nonduplicates(out_file):
    freq_track = defaultdict(int)
    duplicates = []
    nondup = []

    for item in out_file:
        freq_track[item[0] + item[1]] += 1

    for item in out_file:
        if freq_track[item[0] + item[1]] > 1:
            duplicates.append(item)
        elif freq_track[item[0] + item[1]] == 1:
            nondup.append(item)

    return duplicates, nondup

def process_uniq_duplicates(duplicates):
    uniq = []
    sorted_duplicates = sorted(duplicates, key=lambda x: x[1])

    for a, b in zip(sorted_duplicates, sorted_duplicates[1:]):
        if a[1] == b[1]:
            avg_pause = (a[3] + b[3]) / 2
            avg_cov = (a[4] + b[4]) / 2
            uniq.append([a[0], a[1], a[2], avg_pause, avg_cov, a[5], a[6]])

    return uniq

def zscore_calculation(final_uniq):
    lines, zscore_values = [], []
    total_sum = 0
    zscore_hash = {}
    window = 0
    sorted_final_uniq = sorted(final_uniq, key=lambda x: -x[4])

    if len(sorted_final_uniq) >= 300:
        window = 300
    elif len(sorted_final_uniq) == 0:
        print("No pauses predicted for these parameters. Please change the parameters and try again. Thanks.")
        sys.exit()
    else:
        window = len(sorted_final_uniq)

    no_windows = (len(sorted_final_uniq) / window) * 2
    count, i, total, beginat, endat = 1, 1, 0, 1, window

    while count <= no_windows:
        values = []
        total = 0
        for j in range(beginat - 1, min(endat, len(sorted_final_uniq))):
            values.append([sorted_final_uniq[j]])
            total += sorted_final_uniq[j][3]

        count += 1
        mean = total / window
        sqtotal = sum([(mean - item[0][3]) ** 2 for item in values])
        variance = sqtotal / window
        std = math.sqrt(sqtotal / window)

        for k in values:
            for m in k:
                zscore = (m[3] - mean) / std
                zscore_values.append([m[0], m[1], m[2], m[3], m[4], m[5], m[6], zscore])

        beginat = endat - int(window * 0.5)
        endat = beginat + window

    h, uniq_zscore_values = defaultdict(int), {}

    for item in zscore_values:
        if not h[item[0] + item[1]]:
					  h[item[0] + item[1]] += 1
            uniq_zscore_values[item[0] + item[1]] = item
        else:
            uniq_zscore_values[item[0] + item[1]][7] = (uniq_zscore_values[item[0] + item[1]][7] + item[7]) / h[item[0] + item[1]]

    return uniq_zscore_values

def write_output_to_file(outfile, uniq_zscore_values):
    with open(outfile, 'w') as fh:
        for key in sorted(uniq_zscore_values.keys()):
            if key != '':
                fh.write(','.join(str(x) for x in uniq_zscore_values[key]) + '\n')
        print(f"Output has been written to file {outfile}")

def process_annotation(uniq_zscore_values, annotation, db):
    outfile_codon = input("Please enter codon level output file name\n")
    with open(outfile_codon, 'w') as anno:
        anno.write("##gene_name,coordinate_position,number_of_reads_mapped,Pause_score,coverage(%),50_upstream_seq,50_downstream_seq(including_pause_position),Z-score,pause_codon\n")

        with open(annotation, 'r') as f3:
            for line in f3:
                array = line.split()
                annotation_data = {array[0]: [array[1], array[2]]}

                for key in uniq_zscore_values.keys():
                    if uniq_zscore_values[key][0] == array[0]:
                        start, end = int(array[1]), int(array[2])

                        for i in range(start, end + 1, 3):
                            codon_start, codon_end = i, i + 2

                            if codon_start <= uniq_zscore_values[key][1] <= codon_end:
                                seq = db.seq(array[0], codon_start, codon_end)
                                anno.write(','.join(str(x) for x in uniq_zscore_values[key]) + ',' + seq + '\n')


# Main script
if __name__ == "__main__":
    f1_path = input("Enter the F1 file path: ")
    readlength = input("Enter the readlength: ")
    offset = input("Enter the offset: ")
    outfile = input("Enter the output file name: ")
    annotation = input("Enter the annotation file path (leave empty if not available): ")
		fasta_file = input("Enter the FASTA file containing gene sequences: ")
    db = load_sequences(fasta_file)

    with open(f1_path, 'r') as f1:
        final_data = process_file(f1, readlength, offset)

    duplicates, nondup = process_duplicates_nonduplicates(final_data)
    uniq = process_uniq_duplicates(duplicates)
    final_uniq = uniq + nondup
    uniq_zscore_values = zscore_calculation(final_uniq)
    write_output_to_file(outfile, uniq_zscore_values)

    if annotation:
        process_annotation(uniq_zscore_values, annotation, db)
