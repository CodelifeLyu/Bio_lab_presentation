# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 21:19:15 2023

@author: dongs
"""
import random
import math

def generate_dna_strand(length):
    bases = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(bases) for _ in range(length))

def get_complementary_base(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'

def generate_complementary_strand(dna_strand):
    return ''.join(get_complementary_base(base) for base in dna_strand)

def print_dna_structure(strand1, strand2):
    print("5'-", end="")
    for base1 in strand1:
        print(base1, end="")
    print("-3'")

    print("3'-", end="")
    for base2 in strand2:
        print(base2, end="")
    print("-5'")

def reverse_complement(dna_strand):
    complement = generate_complementary_strand(dna_strand)
    return complement[::-1]

def transcribe(dna_strand):
    return dna_strand.replace("T", "U")

def get_amino_acid(codon):
    codon_table = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
        "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    return codon_table.get(codon, "?")

def translate(mrna):
    protein = []
    for i in range(0, len(mrna) - 2, 3):
        codon = mrna[i:i+3]
        amino_acid = get_amino_acid(codon)
        if amino_acid == "*":
            break
        protein.append(amino_acid)
    return ''.join(protein)


length = 9  # Adjust the length of the DNA strand as desired
dna_strand1 = generate_dna_strand(length)
dna_strand2 = generate_complementary_strand(dna_strand1)

print("Original Strand:    ", dna_strand1)
print("Complementary Strand:", dna_strand2)
print("\nDNA Structure:")
print_dna_structure(dna_strand1, dna_strand2)
mRNA_strand = transcribe(dna_strand1)
protein = translate(mRNA_strand)
print("mRNA Strand:              ", mRNA_strand)
print("Protein:                  ", protein)
