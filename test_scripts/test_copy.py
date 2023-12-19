#!/usr/local/bin/python3

"""Copy that contains the print statements"""

from Bio import SeqIO
from math import log
from codon_table import codon_table

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

# Count the number of synonymous substitution and nonsynonymous substitutions between the two sequences
def count_syn_nonsyn_subs(seq1, seq2):
    syn_subs_total = 0
    nonsyn_subs_total = 0
    syn_subs_seq1 = 0
    nonsyn_subs_seq1 = 0
    syn_subs_seq2 = 0
    nonsyn_subs_seq2 = 0

    for codon1, codon2 in zip(split_sequences(seq1), split_sequences(seq2)):
        if is_valid_codon(codon1) and is_valid_codon(codon2):
            aa1 = codon_table.get(codon1, "")
            aa2 = codon_table.get(codon2, "")               
            if codon1 != codon2:
                if aa1 != aa2:
                    nonsyn_subs_total += 1
                    nonsyn_subs_seq1 += 1
                    nonsyn_subs_seq2 += 1
                else:
                    syn_subs_total += 1
                    syn_subs_seq1 += 1
                    syn_subs_seq2 += 1

    print("\nSummary:")
    print(f"Total Synonymous substitutions (seq1): {round(syn_subs_seq1, 3)}")
    print(f"Total Nonsynonymous substitutions (seq1): {round(nonsyn_subs_seq1, 3)}")
    print(f"Total Synonymous substitutions (seq2): {round(syn_subs_seq2, 3)}")
    print(f"Total Nonsynonymous substitutions (seq2): {round(nonsyn_subs_seq2, 3)}")
    print(f"Total Synonymous substitutions: {round(syn_subs_total, 3)}")
    print(f"Total Nonsynonymous substitutions: {round(nonsyn_subs_total, 3)}")

    return syn_subs_seq1, nonsyn_subs_seq1, syn_subs_seq2, nonsyn_subs_seq2, syn_subs_total, nonsyn_subs_total

# Get all possible 1 nucleotide mutant codons at each nucleotide position
def get_possible_mutations(codon, position):
    mutations = [codon[:position] + new_base + codon[position + 1:] for new_base in 'ACGT' if new_base != codon[position]]
    print(mutations)
    return mutations

# Count the number of synonymous and nonsynonymous sites across all codons
def count_syn_nonsyn_sites(seq1, seq2):
    syn_sites_total = 0
    nonsyn_sites_total = 0
    syn_sites_seq1 = 0
    nonsyn_sites_seq1 = 0
    syn_sites_seq2 = 0
    nonsyn_sites_seq2 = 0

    for sites, seq in enumerate([seq1, seq2]):
        print(f"For sequence {sites + 1}:")
        for codon in split_sequences(seq):
            print(f"For codon ({codon}):")
            nonsyn_sites_codon = 0

            for position in range(3):
                mutations = get_possible_mutations(codon, position)
                aa_original = codon_table.get(codon, '')
                aa_mutations = [codon_table.get(mutant_codon, '') for mutant_codon in mutations]

                print(f"  Position {position + 1} ({codon[position]}):")
                for mutant_codon, aa_mutant in zip(mutations, aa_mutations):
                    print(f"    Mutant Codon: {mutant_codon}")
                    print(f"    AA Original: {aa_original}")
                    print(f"    AA Mutant: {aa_mutant} ({'Nonsynonymous' if aa_mutant != aa_original else 'Synonymous'})")

                    if aa_mutant != aa_original:
                        nonsyn_sites_codon += 1 / 3

            syn_sites_codon = 3 - nonsyn_sites_codon
            print(f"  AA Original for {codon}: {aa_original}")
            print(f"Nonsynonymous sites for {codon}: {round(nonsyn_sites_codon,3)}")
            print(f"Synonymous sites for {codon}: {syn_sites_codon}")

            if sites == 0:
                nonsyn_sites_seq1 += nonsyn_sites_codon
                syn_sites_seq1 += syn_sites_codon
            else:
                nonsyn_sites_seq2 += nonsyn_sites_codon
                syn_sites_seq2 += syn_sites_codon

        if sites == 0:
            nonsyn_sites_total += nonsyn_sites_seq1
            syn_sites_total += syn_sites_seq1
        else:
            nonsyn_sites_total += nonsyn_sites_seq2
            syn_sites_total += syn_sites_seq2

    print("\nSummary:")
    print(f"Total Synonymous sites (seq1): {round(syn_sites_seq1,3)}")
    print(f"Total Nonsynonymous sites (seq1): {round(nonsyn_sites_seq1, 3)}")
    print(f"Total Synonymous sites (seq2): {round(syn_sites_seq2, 3)}")
    print(f"Total Nonsynonymous sites (seq2): {round(nonsyn_sites_seq2, 3)}")
    print(f"Total Synonymous sites: {round(syn_sites_total, 3)}")
    print(f"Total Nonsynonymous sites: {round(nonsyn_sites_total, 3)}")
    print(f"Total codons: {len(split_sequences(seq1)) + len(split_sequences(seq2))}")
    
    return (
        syn_sites_seq1, nonsyn_sites_seq1,
        syn_sites_seq2, nonsyn_sites_seq2,
        syn_sites_total, nonsyn_sites_total
    )

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

def calc_pnps(seq1, seq2):
    seq1, seq2 = prep_sequences(seq1, seq2)

    (syn_subs_seq1, nonsyn_subs_seq1, syn_subs_seq2, nonsyn_subs_seq2,
    syn_subs_total, nonsyn_subs_total) = count_syn_nonsyn_subs(seq1, seq2)

    (syn_sites_seq1, nonsyn_sites_seq1, syn_sites_seq2, nonsyn_sites_seq2, 
    syn_sites_total, nonsyn_sites_total) = count_syn_nonsyn_sites(seq1, seq2)

    if syn_sites_total == 0 or nonsyn_sites_total == 0:
        print("Error: Cannot divide by 0, synonymous sites or nonsynonymous sites is 0")
        return 0, 0
    
    pn = (nonsyn_subs_total / nonsyn_sites_total)
    ps = (syn_subs_total / syn_sites_total)

    print("pn:", round(pn, 3))
    print("ps:", round(ps, 3))
    return pn, ps

def calc_dnds_ratio(seq1, seq2):
    pn, ps = calc_pnps(seq1, seq2)

    dn = -(3 / 4) * log(1 - (4 * pn / 3))
    ds = -(3 / 4) * log(1 - (4 * ps / 3))

    if ds == 0 or dn == 0:
        print("Error: dN or dS is 0, cannot calculate dN/dS ratio")
        return 0 
    dnds = dn / ds
    print("dN/dS Ratio:", round(dnds, 3))
    return dnds

aligned_file_path = "/var/www/html/mnguye87/final/test_pwa_files/test_example_2.fasta"
seq1, seq2 = parse_fasta(aligned_file_path)
dnds_ratio = calc_dnds_ratio(seq1, seq2)