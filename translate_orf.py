#! /usr/bin/env python3

"""
Find the first open-reading frame (ORF) in a nucleotide sequence and
translate it into an amino acid sequence.

This is modeled on find_orf.py -- the sequence can be passed directly
on the command line, or read from a file. Start
and stop codons can be customized
"""

import sys

from find_orf import find_first_orf, parse_sequence_from_path
from translate import translate_sequence


# Standard genetic code copied from translate.py

STANDARD_GENETIC_CODE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def main():
    import argparse

    parser = argparse.ArgumentParser(
            description = ('Find the first open-reading frame in a '
                           'nucleotide sequence and print its translation.'))

    default_start_codons = ['AUG']
    default_stop_codons = ['UAA', 'UAG', 'UGA']

    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument should be treated as a path to '
                    'a file containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codon',
            type = str,
            action = 'append',
            default = None,
            help = ('A start codon. This option can be used multiple times '
                    'if there are multiple start codons. '
                    'Default: {0}.'.format(' '.join(default_start_codons))))
    parser.add_argument('-x', '--stop-codon',
            type = str,
            action = 'append',
            default = None,
            help = ('A stop codon. This option can be used multiple times '
                    'if there are multiple stop codons. '
                    'Default: {0}.'.format(' '.join(default_stop_codons))))

    args = parser.parse_args()

    # Pull the sequence from a file or straight from the CLI
    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    # Fall back to default codons if none were supplied
    if not args.start_codon:
        args.start_codon = default_start_codons
    if not args.stop_codon:
        args.stop_codon = default_stop_codons

    # Find the first ORF in the sequence
    orf = find_first_orf(sequence = sequence,
            start_codons = args.start_codon,
            stop_codons = args.stop_codon)

    # Translate the ORF. If no ORF was found, orf is empty string and
    # translate_sequence will return an empty string as well.
    amino_acid_sequence = translate_sequence(
            rna_sequence = orf,
            genetic_code = STANDARD_GENETIC_CODE)

    sys.stdout.write('{}\n'.format(amino_acid_sequence))


if __name__ == '__main__':
    main()
