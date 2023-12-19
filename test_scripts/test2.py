#!/usr/local/bin/python3

"""Checking to see if sequences are being prepepd properly for calculations"""

from Bio import SeqIO

def is_valid_codon(codon):
    valid_bases = {'A', 'T', 'G', 'C'}
    return all(base.upper() in valid_bases for base in ''.join(codon))

def parse_fasta(aligned_file_path):
    sequences = list(SeqIO.parse(aligned_file_path, "fasta"))

    if len(sequences) < 2:
        raise ValueError("For pairwise comparison, at least two sequences are required for dN/dS calculation.")
    
    seq1 = str(sequences[0].seq)
    seq2 = str(sequences[1].seq)
    return seq1, seq2

def clean_sequences(seq1, seq2):
    seq1 = ''.join(c.replace('-', '').replace(' ', '') for c in seq1)
    seq2 = ''.join(c.replace('-', '').replace(' ', '') for c in seq2)  
    return seq1, seq2

def trim_sequences(seq):
    if len(seq) % 3 != 0:
        adjustment = len(seq) % 3
        print(f"Adjustment: Removed {adjustment} base(s)")
        seq = seq[:-adjustment]
        print("Adjusted Length:", len(seq))
    return seq

def split_sequences(seq):
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]

def prep_sequences(seq1, seq2):
    seq1, seq2 = clean_sequences(seq1, seq2)
    print("Cleaned Length of seq1:", len(seq1))
    print("Cleaned Length of seq2:", len(seq2))

    if len(seq1) != len(seq2):
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]
        print(f"Trimmed to equal length: {min_len}")

    seq1 = trim_sequences(seq1)
    seq2 = trim_sequences(seq2)
        
    print("Trimmed Length of seq1:", len(seq1))
    print("Trimmed Length of seq2:", len(seq2))
    return seq1, seq2

aligned_file_path = "/var/www/html/mnguye87/final/test_pwa_files/OPN1LW_cat_human.fasta"
seq1, seq2 = parse_fasta(aligned_file_path)

print("Original Length of seq1:", len(seq1))
print("Original Length of seq2:", len(seq2))

seq1, seq2 = prep_sequences(seq1, seq2)

print("Processed Length of seq1:", len(seq1))
print("Processed Length of seq2:", len(seq2))

# Summary:
# Total Synonymous sites: 6943.3333333330365
# Total Nonsynonymous sites: 26350.66666666764
# Total codons: 11098
# pn: 0.198
# ps: 0.033
# dN/dS Ratio: 6.835